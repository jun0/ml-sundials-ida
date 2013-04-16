(* Low-level interface to IDA.
   TODO: package this into a less error-prone interface which, perhaps, supplies
   all the data (solver, tolerance, etc.) in one go instead of in pieces.
 *)

(* Static configuration of SUNDIALS.  *)
module SUNDIALS =
struct
  open Bigarray
  type realtype = float64_elt
end

(* IDA solver with serial nvectors and dense matrices.  TODO: abstract out the
   types of vectors and of matrices, so that different configurations can be
   built up by functor application.  (But to do that I need to first get to
   know the range of configurations possible in C.)  *)
module SerialDense =
struct
  open Bigarray
  include SUNDIALS

  (* Serial N_Vector is represented as a bigarray.  *)
  type nvector = (float, realtype, c_layout) Array1.t
  let nvector_of_array = Array1.of_array float64 c_layout
  let nvector_of_list ls = nvector_of_array (Array.of_list ls)
  (* Dense DlsMat is represented as a bigarray.  *)
  type dlsmat = (float, realtype, c_layout) Array2.t

  type residual_fn = float -> nvector -> nvector -> nvector -> int
  type ida
  external ida_init : residual_fn -> float -> nvector -> nvector -> ida =
    "caml_ida_init"
  external ida_sv_tolerances : ida -> float -> nvector -> unit =
    "caml_ida_sv_tolerances"

  (* TODO: allow the user to throw an exception instead of returning
     JacFnError.  *)
  type ida_dls_jac_fn_result = JacFnOK | JacFnStepTooBig | JacFnError
  type ida_dls_jac_fn =
    float ->                            (* tt *)
    float ->                            (* cj *)
    nvector ->                          (* yy *)
    nvector ->                          (* yp (i.e. y'(t)) *)
    nvector ->                          (* rr *)
    dlsmat ->                           (* jacobian (output) *)
    nvector ->                          (* scratch space *)
    nvector ->                          (* scratch space *)
    nvector ->                          (* scratch space *)
    ida_dls_jac_fn_result

  (* Attach a direct linear solver with the given Jacobian function.  *)
  external ida_jacobian_fn : ida -> ida_dls_jac_fn -> unit
    = "caml_ida_dense_jacobian_fn"

  external ida_get_last_order : ida -> int = "caml_ida_get_last_order"
  external ida_get_num_steps : ida -> int = "caml_ida_get_num_steps"
  external ida_get_last_step : ida -> float = "caml_ida_get_last_step"
  external ida_get_num_res_evals : ida -> int = "caml_ida_get_num_res_evals"
  external ida_dls_get_num_res_evals : ida -> int =
    "caml_ida_dls_get_num_res_evals"
  external ida_dls_get_num_jac_evals : ida -> int =
    "caml_ida_dls_get_num_jac_evals"
  external ida_get_num_nonlin_solv_iters : ida -> int =
    "caml_ida_get_num_nonlin_solv_iters"
  external ida_get_num_nonlin_solv_conv_fails : ida -> int =
    "caml_ida_get_num_nonlin_solv_conv_fails"
  external ida_get_num_err_test_fails : ida -> int =
    "caml_ida_get_num_err_test_fails"
  external ida_get_num_g_evals : ida -> int = "caml_ida_get_num_g_evals"

  external ida_get_root_info : ida -> int array -> unit =
    "caml_ida_get_root_info"

  type solve_result = IDA_SUCCESS | IDA_ROOT_RETURN of float
  type solve_task = IDA_NORMAL | IDA_ONE_STEP
  external ida_solve : ida -> float -> nvector -> nvector ->
    solve_task -> solve_result = "caml_ida_solve"

  type root_fn =
    float ->                            (* t *)
    nvector ->                          (* y *)
    nvector ->                          (* y' *)
    nvector ->                          (* put values here *)
    int                                 (* 0 for success nonzero for failure *)
  external ida_root_init : ida -> int -> root_fn -> unit = "caml_ida_root_init"

end
