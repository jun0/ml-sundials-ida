#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/custom.h>
#include <caml/bigarray.h> 
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <ida/ida_dense.h>
#include <sundials/sundials_types.h>

#if defined (SUNDIALS_DOUBLE_PRECISION)
#define CAML_SUNDIALS_BA_KIND CAML_BA_FLOAT64
#else
#error "ML-SUNDIALS requires SUNDIALS to be compiled with double precision."
#endif

#include <assert.h>
#include <stdlib.h>

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}


struct ida {
  void *mem;
  int num_eqs;
  int num_root_fns;
  int *root_info_buf;
  value caml_data;
};

enum ida_caml_data_indices {
  CAML_DATA_RESIDUAL_FN = 0,
  CAML_DATA_JAC_FN,
  CAML_DATA_ROOT_FN,
  CAML_DATA_RESIDUAL_TRES,
  CAML_DATA_RESIDUAL_YY,
  CAML_DATA_RESIDUAL_YP,
  CAML_DATA_RESIDUAL_RR,

  CAML_DATA_RESIDUAL_DATA_SIZE
};

/* Access a struct ida pointer (as lvalue) stored in a custom OCaml value.  */
#define IDA_val(v) (*(struct ida**) Data_custom_val (v))

/* Wraps a bigarray's storage in a serial N_Vector.  The resulting N_Vector
 * shares storage with the bigarray, but does not take ownership of it.  The
 * N_Vector must be explicitly destroyed with N_VDestroy_Serial() before the
 * bigarray has a chance to be collected by OCaml's GC; the N_VDestroy_Serial()
 * will not deallocate any portion of the bigarray.  */
N_Vector nvector_of_bigarray (value v)
{
  assert (Bigarray_val (v)->num_dims == 1);
  return N_VMake_Serial (Bigarray_val(v)->dim[0],
                         Data_bigarray_val (v));
}

/* Wraps a serial N_Vector's storage in a bigarray.  The resulting bigarray
 * shrares storage with the N_Vector but does not take ownership of it.  If the
 * original N_Vector claimed ownership of the storage, then the bigarray must
 * not be accessed after the N_Vector is destroyed with
 * N_VDestroy_Serial().  */
value bigarray_of_nvector (N_Vector v)
{
  CAMLparam0 ();
  CAMLlocal1 (ret);
  intnat dims[1] = { NV_LENGTH_S (v) };
  ret = caml_ba_alloc (CAML_BA_C_LAYOUT | CAML_SUNDIALS_BA_KIND
                       | CAML_BA_EXTERNAL,
                       1, NV_DATA_S (v), dims);
  CAMLreturn (ret);
}

/* Wraps a C array in a bigarray.  The resulting bigarray does not take
 * ownership of the C array.  The bigarray must not be accessed after the C
 * array has been deallocated.  */
value bigarray_of_array (realtype *buf, long length)
{
  CAMLparam0 ();
  CAMLlocal1 (ret);
  intnat dims[1] = { length };
  ret = caml_ba_alloc (CAML_BA_C_LAYOUT | CAML_SUNDIALS_BA_KIND
                       | CAML_BA_EXTERNAL,
                       1, buf, dims);
  CAMLreturn (ret);
}

/* Wraps a dense DlsMat's storage in a bigarray.  The resulting bigarray
 * shrares storage with the DlsMat but does not take ownership of it.  The
 * bigarray must not be accessed after the DlsMat is destroyed.  */
value bigarray_of_dense_dlsmat (DlsMat m)
{
  CAMLparam0 ();
  CAMLlocal1 (ret);
  intnat dims[2] = { m->N, m->M }; /* num_cols, num_rows */
  ret = caml_ba_alloc (CAML_BA_C_LAYOUT | CAML_SUNDIALS_BA_KIND
                       | CAML_BA_EXTERNAL,
                       2, m->data, dims);
  CAMLreturn (ret);
}


void caml_destroy_ida (value val_ida)
{
  struct ida *ida = IDA_val (val_ida);
  IDAFree (&ida->mem);
  caml_remove_global_root (&ida->caml_data);
  free (ida->root_info_buf);
  free (ida);
}

static struct custom_operations caml_ida_ops = {
  "ida",
  caml_destroy_ida,
  custom_compare_default,
  custom_hash_default,
  custom_serialize_default,
  custom_deserialize_default

  /* FIXME: this field is present only for OCaml >= 3.12.1.  How do I check
   * OCaml's version?  */
  , custom_compare_default        /* compare_ext; not explained in manual */
};

int callback_residual_fn (realtype tres, N_Vector yy, N_Vector yp,
                          N_Vector rr, void *user_data)
{
  CAMLparam0 ();
  struct ida *ida = user_data;
  CAMLlocalN (args, 4);
  intnat dims[1] = {0};

  /* Convert everything to OCaml values.  */
  args[0] = Field (ida->caml_data, CAML_DATA_RESIDUAL_TRES);
  Field (args[0], 0) = tres;
  dims[0] = NV_LENGTH_S (yy);
  args[1] = caml_ba_alloc (CAML_BA_C_LAYOUT | CAML_SUNDIALS_BA_KIND
                           | CAML_BA_EXTERNAL,
                           1, NV_DATA_S (yy), dims);
  dims[0] = NV_LENGTH_S (yp);
  args[2] = caml_ba_alloc (CAML_BA_C_LAYOUT | CAML_SUNDIALS_BA_KIND
                           | CAML_BA_EXTERNAL,
                           1, NV_DATA_S (yp), dims);
  dims[0] = NV_LENGTH_S (rr);
  args[3] = caml_ba_alloc (CAML_BA_C_LAYOUT | CAML_SUNDIALS_BA_KIND
                           | CAML_BA_EXTERNAL,
                           1, NV_DATA_S (rr), dims);

  caml_callbackN (Field (ida->caml_data, CAML_DATA_RESIDUAL_FN), 4, args);
  CAMLreturnT (int, 0);
}

/* ida_init : residual_fn -> float -> vector -> vector -> ida */
CAMLprim value caml_ida_init (value val_residual_fn, value val_t0,
                              value val_yy, value val_yp)
{
  CAMLparam4 (val_residual_fn, val_t0, val_yy, val_yp);
  CAMLlocal2 (caml_data_tuple, ida);
  N_Vector yy = nvector_of_bigarray (val_yy);
  N_Vector yp = nvector_of_bigarray (val_yp);
  double t0 = Double_val (val_t0);
  int dims[] = {0};
  int retflag;

  struct ida *ret = (struct ida*) malloc (sizeof (struct ida));
  ret->mem = IDACreate ();
  ret->num_eqs = NV_LENGTH_S (yy);
  ret->num_root_fns = 0;
  ret->root_info_buf = NULL;
  /* TODO: fail if yy and yp differ in length.  */

  if (check_flag (ret->mem, "IDACreate", 0))
    abort ();                   /* TODO: raise exception */
  /* TODO: error check */
  retflag = IDAInit (ret->mem, callback_residual_fn, t0, yy, yp);
  if (check_flag (&retflag, "IDAInit", 1))
    abort ();                   /* TODO: raise exception */

  caml_register_global_root (&ret->caml_data);
  ret->caml_data = caml_alloc_tuple (CAML_DATA_RESIDUAL_DATA_SIZE);

  Store_field (ret->caml_data, CAML_DATA_RESIDUAL_FN, val_residual_fn);
  Store_field (ret->caml_data, CAML_DATA_RESIDUAL_TRES, caml_copy_double (0.0));
  /* TODO: avoid allocation of bigarrays during callback.  */

  retflag = IDASetUserData (ret->mem, ret);
  if (check_flag (&retflag, "IDASetUserData", 1))
    abort ();                   /* TODO: raise exception */

  /* Wrap ret in an OCaml object.  Try to have no more than ~10M vector
   * elements worth of floating garbage.  */
  ida = caml_alloc_custom (&caml_ida_ops,
                           sizeof (ret),
                           Bigarray_val (val_yy)->dim[0],
                           10 * 1024 * 1024);
  IDA_val (ida) = ret;

  /* IDAInit copies everything it needs from the vectors it receives, so we're
   * done with yy and yp.  */
  N_VDestroy_Serial (yy);
  N_VDestroy_Serial (yp);

  CAMLreturn (ida);
}

/* ida_sv_tolerances : ida -> float -> vector -> unit */
CAMLprim void caml_ida_sv_tolerances (value val_ida, value val_rtol,
                                      value val_avtol)
{
  CAMLparam3 (val_ida, val_rtol, val_avtol);
  struct ida *ida = IDA_val (val_ida);
  double rtol = Double_val (val_rtol);
  N_Vector avtol = nvector_of_bigarray (val_avtol);
  int retval;

  retval = IDASVtolerances (ida->mem, rtol, avtol);
  if (check_flag (&retval, "IDASVtolerances", 1))
    abort ();                   /* TODO: raise exception */

  N_VDestroy_Serial (avtol);
  CAMLreturn0;
}

int callback_dense_jacobian_fn (long int Neq, realtype tt,  realtype cj, 
                                N_Vector yy, N_Vector yp, N_Vector resvec,
                                DlsMat JJ, void *user_data,
                                N_Vector tempv1, N_Vector tempv2,
                                N_Vector tempv3)
{
#define NARGS 9
  CAMLparam0 ();
  struct ida *ida = user_data;
  CAMLlocalN (args, NARGS);
  CAMLlocal1 (retcode);
  intnat dims[1] = {0};

  /* Convert everything to OCaml values.  */
  args[0] = caml_copy_double (tt);
  args[1] = caml_copy_double (cj);
  args[2] = bigarray_of_nvector (yy);
  args[3] = bigarray_of_nvector (yp);
  args[4] = bigarray_of_nvector (resvec);
  args[5] = bigarray_of_dense_dlsmat (JJ);
  args[6] = bigarray_of_nvector (tempv1);
  args[7] = bigarray_of_nvector (tempv2);
  args[8] = bigarray_of_nvector (tempv3);

  retcode = caml_callbackN (Field (ida->caml_data, CAML_DATA_JAC_FN),
                            NARGS, args);
#define JAC_FN_OK (Val_int (0))
#define JAC_FN_STEP_TOO_BIG (Val_int (1))
#define JAC_FN_ERROR (Val_int (2))
  switch (retcode)
    {
    case JAC_FN_OK:
      CAMLreturnT (int, 0);
    case JAC_FN_STEP_TOO_BIG:
      CAMLreturnT (int, 1);
    default:
      CAMLreturnT (int, -1);
    }
#undef NARGS
}


/* caml_ida_dense_jacobian_fn : ida -> ida_dls_jac_fcn -> int */
CAMLprim void caml_ida_dense_jacobian_fn (value val_ida, value val_jac_fn)
{
  CAMLparam2 (val_ida, val_jac_fn);
  int retval;
  struct ida *ida = IDA_val (val_ida);

  retval = IDADense (ida->mem, ida->num_eqs);
  if (check_flag (&retval, "IDADense", 1))
    abort ();                   /* TODO: raise an exception */

  retval = IDADlsSetDenseJacFn (ida->mem, callback_dense_jacobian_fn);
  if (check_flag (&retval, "IDADlsSetDenseJacFn", 1))
    abort ();                   /* TODO: raise an exception */

  Store_field (ida->caml_data, CAML_DATA_JAC_FN, val_jac_fn);

  CAMLreturn0;
}


int callback_root_fn (realtype tt, N_Vector yy, N_Vector yp, realtype *gout,
                      void *user_data)
{
#define NARGS 4
  CAMLparam0 ();
  struct ida *ida = user_data;
  CAMLlocalN (args, NARGS);
  CAMLlocal1 (retcode);
  intnat dims[1] = {0};

  /* Convert everything to OCaml values.  */
  args[0] = caml_copy_double (tt);
  args[1] = bigarray_of_nvector (yy);
  args[2] = bigarray_of_nvector (yp);
  args[3] = bigarray_of_array (gout, ida->num_root_fns);

  retcode = caml_callbackN (Field (ida->caml_data, CAML_DATA_ROOT_FN),
                            NARGS, args);
  CAMLreturnT (int, Int_val (retcode));
}


/* caml_ida_root_init : ida -> root_fcn -> int -> unit */
CAMLprim void caml_ida_root_init (value val_ida, value val_nrtfn,
                                  value val_root_fn)
{
  CAMLparam3 (val_ida, val_nrtfn, val_root_fn);
  int retval;
  struct ida *ida = IDA_val (val_ida);
  int nrtfn = Int_val (val_nrtfn);

  retval = IDARootInit (ida->mem, nrtfn, callback_root_fn);
  if (check_flag (&retval, "IDARootInit", 1))
    abort ();                   /* TODO: raise an exception */

  Store_field (ida->caml_data, CAML_DATA_ROOT_FN, val_root_fn);
  ida->num_root_fns = nrtfn;
  ida->root_info_buf = realloc (ida->root_info_buf,
                                sizeof (ida->root_info_buf[0]) * nrtfn);

  CAMLreturn0;
}


#define CAML_IDA_NORMAL (Val_int(0))
#define CAML_IDA_ONE_STEP (Val_int(1))

#define CAML_IDA_SUCCESS (Val_int(0))
#define CAML_IDA_ROOT_RETURN_TAG 0
value caml_ida_root_return (double tret)
{
  CAMLparam0 ();
  CAMLlocal2 (ret, f);
  f = caml_copy_double (tret);
  ret = caml_alloc (1, CAML_IDA_ROOT_RETURN_TAG);
  Store_field (ret, 0, f);
  CAMLreturn (ret);
}

int task_tag_to_int (value task)
{
  assert (task == CAML_IDA_NORMAL || task == CAML_IDA_ONE_STEP);
  if (task == CAML_IDA_NORMAL)
    return IDA_NORMAL;
  return IDA_ONE_STEP;
}

CAMLprim value caml_ida_solve (value val_ida, value val_tout,
                               value val_yy, value val_yp, value task)
{
  CAMLparam5 (val_ida, val_tout, val_yy, val_yp, task);
  CAMLlocal1 (ret);
  int retcode;
  struct ida *ida = IDA_val (val_ida);
  double tout = Double_val (val_tout);
  N_Vector yy = nvector_of_bigarray (val_yy);
  N_Vector yp = nvector_of_bigarray (val_yp);
  double tret;

  retcode = IDASolve (ida->mem, tout, &tret, yy, yp, task_tag_to_int (task));

  if (check_flag(&retcode, "IDASolve", 1))
    abort ();
    
  switch (retcode)
    {
    case IDA_SUCCESS:
      ret = CAML_IDA_SUCCESS;
      break;
    case IDA_ROOT_RETURN:
      ret = caml_ida_root_return (tret);
      break;
    default:
      assert (FALSE);
    }

  N_VDestroy_Serial (yy);
  N_VDestroy_Serial (yp);
  
  CAMLreturn (ret);
}


/* caml_ida_get_last_order : ida -> int */
CAMLprim value caml_ida_get_last_order (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  int ret;
  IDAGetLastOrder (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_int (ret));
}

/* caml_ida_get_num_steps : ida -> int */
CAMLprim value caml_ida_get_num_steps (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDAGetNumSteps (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_last_step : ida -> float */
CAMLprim value caml_ida_get_last_step (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  double ret;
  IDAGetLastStep (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (caml_copy_double (ret));
}

/* caml_ida_dls_get_num_res_evals : ida -> int */
CAMLprim value caml_ida_dls_get_num_res_evals (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDADlsGetNumResEvals (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_num_res_evals : ida -> int */
CAMLprim value caml_ida_get_num_res_evals (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDAGetNumResEvals (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_dls_get_num_jac_evals : ida -> int */
CAMLprim value caml_ida_dls_get_num_jac_evals (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDADlsGetNumJacEvals (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_num_nonlin_solv_iters : ida -> int */
CAMLprim value caml_ida_get_num_nonlin_solv_iters (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDAGetNumNonlinSolvIters (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_num_nonlin_solv_conv_fails : ida -> int */
CAMLprim value caml_ida_get_num_nonlin_solv_conv_fails (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDAGetNumNonlinSolvConvFails (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_num_err_test_fails : ida -> int */
CAMLprim value caml_ida_get_num_err_test_fails (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDAGetNumErrTestFails (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_num_g_evals : ida -> int */
CAMLprim value caml_ida_get_num_g_evals (value val_ida)
{
  CAMLparam1 (val_ida);
  struct ida *ida = IDA_val (val_ida);
  long ret;
  IDAGetNumGEvals (ida->mem, &ret);
  /* TODO: check return value */
  CAMLreturn (Val_long (ret));
}

/* caml_ida_get_root_info : ida -> int array -> unit */
CAMLprim void caml_ida_get_root_info (value val_ida, value val_buf)
{
  CAMLparam2 (val_ida, val_buf);
  struct ida *ida = IDA_val (val_ida);
  int i;

  IDAGetRootInfo (ida->mem, ida->root_info_buf);
  fflush (stdout);
  for (i = 0; i < ida->num_root_fns; ++i)
    Store_field (val_buf, i, Val_int (ida->root_info_buf[i]));

  CAMLreturn0;
}
