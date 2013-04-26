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
  module NVector = struct
    type t = (float, realtype, c_layout) Array1.t
    let of_array = Array1.of_array float64 c_layout
    let of_list ls = of_array (Array.of_list ls)
    let create = Array1.create float64 c_layout
    let blit = Array1.blit
    let fill = Array1.fill
    let get = Array1.get
    let set = Array1.set

    let maxnorm u =
      let rec go x i =
        if i < Array1.dim u then
          go (max x (abs_float u.{i})) (i+1)
        else x
      in go u.{0} 1
  end
  type nvector = NVector.t
  (* Dense DlsMat is represented as a bigarray.  *)
  type dlsmat = (float, realtype, c_layout) Array2.t

  type residual_fn = float -> nvector -> nvector -> nvector -> int
  type ida
  external ida_init : residual_fn -> float -> nvector -> nvector -> ida =
    "caml_ida_init"
  external ida_sv_tolerances : ida -> float -> nvector -> unit =
    "caml_ida_sv_tolerances"
  external ida_ss_tolerances : ida -> float -> float -> unit =
    "caml_ida_ss_tolerances"

  external ida_band : ida -> int -> int -> int -> unit =
    "caml_ida_band"

  type calc_ic_opt = IDA_YA_YDP_INIT | IDA_Y_INIT
  external ida_calc_ic : ida -> calc_ic_opt -> float -> unit =
    "caml_ida_calc_ic"

  external ida_ss_tolerances : ida -> float -> float -> unit =
    "caml_ida_ss_tolerances"

  (* TODO: allow the user to throw an exception instead of returning
     JacFnError.  *)
  type ida_dls_jac_fn_result = JacFnOK | JacFnStepTooBig | JacFnError
  type ida_dls_jac_fn =
    float ->                            (* t *)
    float ->                            (* cj *)
    nvector ->                          (* y(t) *)
    nvector ->                          (* y'(t) *)
    nvector ->                          (* rr *)
    dlsmat ->                           (* jacobian (output) *)
    nvector ->                          (* scratch space *)
    nvector ->                          (* scratch space *)
    nvector ->                          (* scratch space *)
    ida_dls_jac_fn_result

  (* Attach a direct linear solver with the given Jacobian function.  *)
  external ida_jacobian_fn : ida -> ida_dls_jac_fn -> unit
    = "caml_ida_dense_jacobian_fn"

  (* Set algebraic/differential components in the y vector.  Required for
     ida_calc_ic.  *)
  external ida_set_id : ida -> nvector -> unit =
    "caml_ida_set_id"

  (* Specify inequality constraints for solution vector y.  At each dimension,
     the nvector should be:
        0.0 for no constraint
        1.0 for y.{i} >= 0.0
       -1.0 for y.{i} <= 0.0
        2.0 for y.{i} > 0.0
       -2.0 for y.{i} < 0.0
     TODO: provide interface with this data type:
     type ida_constraint =
       UNCONSTRAINED | NON_NEGATIVE | NON_POSITIVE | POSITIVE | NEGATIVE
   *)
  external ida_set_constraints : ida -> nvector -> unit =
    "caml_ida_set_constraints"

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
  external ida_get_num_nonlin_solve_conv_fails : ida -> int =
    "caml_ida_get_num_nonlin_solve_conv_fails"

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
