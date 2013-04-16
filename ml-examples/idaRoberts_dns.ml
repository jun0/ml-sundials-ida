open IDA_Raw.SerialDense
let print_header rtol avtol yy =
  let open Printf in
  printf "\nidaRoberts_dns: Robertson kinetics DAE serial example problem for IDA\n";
  printf "         Three equation chemical kinetics problem.\n\n";
  printf "Linear solver: IDADENSE, with user-supplied Jacobian.\n";

  printf "Tolerance parameters:  rtol = %g   atol = %g %g %g \n"
         rtol avtol.{0} avtol.{1} avtol.{2};
  printf "Initial conditions y0 = (%g %g %g)\n"
         yy.{0} yy.{1} yy.{2};

  printf "Constraints and id not used.\n\n";
  printf "-----------------------------------------------------------------------\n";
  printf "  t             y1           y2           y3";
  printf "      | nst  k      h\n";
  printf "-----------------------------------------------------------------------\n";
and print_output ida t y =
  let open Printf in
  let kused = ida_get_last_order ida
  and nst   = ida_get_num_steps ida
  and hused = ida_get_last_step ida
  in
  printf "%10.4e %12.4e %12.4e %12.4e | %3d  %1d %12.4e\n" 
         t y.{0} y.{1} y.{2} nst kused hused;
and print_final_stats ida =
  let nst = ida_get_num_steps ida
  and nre = ida_get_num_res_evals ida
  and nje = ida_dls_get_num_jac_evals ida
  and nni = ida_get_num_nonlin_solv_iters ida
  and netf = ida_get_num_err_test_fails ida
  and ncfn = ida_get_num_nonlin_solv_conv_fails ida
  and nreLS = ida_dls_get_num_res_evals ida
  and nge = ida_get_num_g_evals ida
  in
  let open Printf in
  printf "\nFinal Run Statistics: \n\n";
  printf "Number of steps                    = %d\n" nst;
  printf "Number of residual evaluations     = %d\n" (nre+nreLS);
  printf "Number of Jacobian evaluations     = %d\n" nje;
  printf "Number of nonlinear iterations     = %d\n" nni;
  printf "Number of error test failures      = %d\n" netf;
  printf "Number of nonlinear conv. failures = %d\n" ncfn;
  printf "Number of root fn. evaluations     = %d\n" nge;
and print_root_info root_f1 root_f2 =
  Printf.printf "    rootsfound[] = %3d %3d\n" root_f1 root_f2;
;;

let resrob tres y yp rr =
  rr.{0} <- -.0.04*.y.{0} +. 1.0e4*.y.{1}*.y.{2};
  rr.{1} <- -.rr.{0} -. 3.0e7*.y.{1}*.y.{1} -. yp.{1};
  rr.{0} <-  rr.{0} -. yp.{0};
  rr.{2} <-  y.{0} +. y.{1} +. y.{2} -. 1.0;
  0
and jacrob tt cj y yp resvec jj _ _ _ =
  (* Note column-major storage.  *)
  jj.{0, 0} <- -.0.04 -. cj;
  jj.{0, 1} <- 0.04;
  jj.{0, 2} <- 1.0;
  jj.{1, 0} <- 1.0e4*.y.{2};
  jj.{1, 1} <- -.1.0e4*.y.{2} -. 6.0e7*.y.{1} -. cj;
  jj.{1, 2} <- 1.0;
  jj.{2, 0} <- 1.0e4*.y.{1};
  jj.{2, 1} <- -.1.0e4*.y.{1};
  jj.{2, 2} <- 1.0;
  JacFnOK
and grob t y y' gout =
  let y1 = y.{0}
  and y3 = y.{2}
  in
  gout.{0} <- y1 -. 0.0001;
  gout.{1} <- y3 -. 0.01;
  0
;;

let y = nvector_of_list [1.; 0.; 0.]
and y' = nvector_of_list [-0.04; 0.04; 0.]
and rtol = 1.0e-4
and avtol = nvector_of_list [1.0e-8; 1.0e-14; 1.0e-6]
and t0 = 0.
and tout1 = 0.4
and nout = 12
and rootsfound = [|0;0|]
in
let _ = print_header rtol avtol y in
let ida = ida_init resrob t0 y y' in
let _ = ida_sv_tolerances ida rtol avtol in
let _ = ida_jacobian_fn ida jacrob in
let _ = ida_root_init ida 2 grob in
let rec loop iout tout =
  Gc.full_major ();
  if iout < nout then
    begin
      match ida_solve ida tout y y' IDA_NORMAL with
      | IDA_SUCCESS -> (print_output ida tout y;
                        loop (iout+1) (tout *. 10.0))
      | IDA_ROOT_RETURN tret ->
        print_output ida tret y;
        ida_get_root_info ida rootsfound;
        print_root_info rootsfound.(0) rootsfound.(1);
        loop iout tout
    end
in
loop 0 tout1;
print_final_stats ida
;;
