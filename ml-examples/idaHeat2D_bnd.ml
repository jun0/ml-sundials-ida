open Bigarray
open Printf
open IDA_Raw.SerialDense

(* Problem Constants *)
let nout = 11
and mgrid = 10
and bval = 0.1
let neq = mgrid * mgrid
;;

let print_header rtol atol =
  printf "\nidaHeat2D_bnd: Heat equation, serial example problem for IDA\n";
  printf "          Discretized heat equation on 2D unit square.\n";
  printf "          Zero boundary conditions,";
  printf " polynomial initial conditions.\n";
  printf "          Mesh dimensions: %d x %d" mgrid mgrid;
  printf "        Total system size: %d\n\n" neq;
  printf "Tolerance parameters:  rtol = %g   atol = %g\n" rtol atol;
  printf "Constraints set to force all solution components >= 0. \n";
  printf "Linear solver: IDABAND, banded direct solver \n";
  printf "       difference quotient Jacobian, half-bandwidths = %d \n" mgrid;
  printf "IDACalcIC called with input boundary values = %g \n" bval;
  (* Print output table heading and initial line of table. *)
  printf "\n   Output Summary (umax = max-norm of solution) \n\n";
  printf "  time       umax     k  nst  nni  nje   nre   nreLS    h      \n" ;
  printf " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n"
and print_output ida t uu =
  let umax = NVector.maxnorm uu
  and kused = ida_get_last_order ida
  and nst = ida_get_num_steps ida
  and nni = ida_get_num_nonlin_solv_iters ida
  and nre = ida_get_num_res_evals ida
  and hused = ida_get_last_step ida
  and nje = ida_dls_get_num_jac_evals ida
  and nreLS = ida_dls_get_num_res_evals ida
  in
  printf " %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e \n"
         t umax kused nst nni nje nre nreLS hused
;;

let print_matrix a m prefix =
  for i = 0 to m-1 do
    printf "%s[" prefix;
    for j = 0 to m-2 do
      printf "%g, " a.{j*m + i}
    done;
    printf "%g]\n" a.{(m-1)*m + i};
  done

let heatres (mm,dx,coeff) tres uu up resval =

  (* Initialize resval to uu, to take care of boundary equations. *)
  Array1.blit uu resval;
  (* Loop over interior points; set res = up - (central difference). *)
  for j = 1 to mm-2 do
    let offset = mm*j in
    for i = 1 to mm-2 do
      let loc = offset + i in
      resval.{loc} <- up.{loc} -. coeff *.
        (uu.{loc-1} +. uu.{loc+1} +. uu.{loc-mm} +. uu.{loc+mm} -. 4.*.uu.{loc})
    done
  done;
  0
let set_initial_profile ((mm,dx,coeff) as params) uu up id res =
  let mm1 = mm - 1 in

  Array1.fill id 1.;
  (* Initialize uu on all grid points.  *)
  for j = 0 to mm1 do
    let yfact = dx *. float_of_int j
    and offset = mm*j
    in
    for i = 0 to mm1 do
      let xfact = dx *. float_of_int i
      and loc = offset + i
      in
      uu.{loc} <- 16. *. xfact *. (1. -. xfact) *. yfact *. (1. -. yfact)
    done
  done;

  Array1.fill up 0.;
  (* heatres sets res to negative of ODE RHS values at interior points.  *)
  ignore (heatres params 0. uu up res);
  (* Copy -res into up to get correct interior initial up values. *)
  for i = 0 to Array1.dim res - 1 do
    up.{i} <- -. res.{i}
  done;

  (* Finally, set values of u, up, and id at boundary points.  *)
  for j = 0 to mm1 do
    let offset = mm*j in
    for i = 0 to mm1 do
      let loc = offset + i in
      if (j = 0 || j = mm1 || i = 0 || i = mm1) then
        (uu.{loc} <- bval;
         up.{loc} <- 0.;
         id.{loc} <- 0.);
    done
  done;

;;

let uu = NVector.create neq
and up = NVector.create neq
and res = NVector.create neq
and constraints = NVector.create neq
and id = NVector.create neq
and dx = 1. /. float_of_int (mgrid - 1)
in
let data = (mgrid, dx, 1. /. (dx *. dx))
and t0   = 0.;
and t1   = 0.01;
and rtol = 0.;
and atol = 1.0e-3;
and mu = mgrid
and ml = mgrid
in
set_initial_profile data uu up id res;

(* Set constraints to all 1's for nonnegative solution values. *)
NVector.fill constraints 1.;

let ida = ida_init (heatres data) t0 uu up in
ida_set_id ida id;
ida_set_constraints ida constraints;
ida_ss_tolerances ida rtol atol;
ida_band ida neq mu ml;
ida_calc_ic ida IDA_YA_YDP_INIT t1;
print_header rtol atol;
print_output ida t0 uu;
let rec go tout iout =
  if iout <= nout then
    begin
      let _ = ida_solve ida tout uu up IDA_NORMAL in
      (* Note: In IDA_NORMAL, tret = tout in the C verison unless there's an
         error, so simply saying `tout' is a faithful translation.  *)
      print_output ida tout uu;
      go (tout *. 2.) (iout + 1)
    end;
in
go t1 1;
let netf = ida_get_num_err_test_fails ida
and ncfn = ida_get_num_nonlin_solv_conv_fails ida
in
printf "\n netf = %d,   ncfn = %d \n" netf ncfn
;Gc.full_major ()
