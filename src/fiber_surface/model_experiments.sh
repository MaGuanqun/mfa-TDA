# Shared-edge queries and resolution settings for the full-domain model experiments.
# The MFA adapter reads both grid extents from the files. CoordNet uses its
# physical input bounds for normalization. Both extract over the resulting domain.
isabel_pressure_temperature() {
    model_args=(--first "${FIRST_MFA:-$fiber_repo/build/src/hurricane_isabel_3d/Pf30_3d.mfa}"
                --second "${SECOND_MFA:-$fiber_repo/build/src/hurricane_isabel_3d/TCf30_3d.mfa}"
                --variable1 0 --variable2 0)
    mesh_args=(--step 1 --spatial-tolerance 0.08 --residual-tolerance 0.005
               --seed-grid 7 --max-triangles 8000)
    edge_a=-1100,-27
    edge_b=-1100,-22
}

vortex_velocities() {
    local model_dir="$fiber_repo/CoordNet/Exp/vortex_street_3d_2_variables"
    model_args=(--first "$model_dir/u-super-spatial-temporal-64-5-float64.pt"
                --second "$model_dir/v-super-spatial-temporal-64-5-float64.pt"
                --domain-min -0.5,-0.5,13.5 --domain-max 7.5,0.5,15)
    vortex_residual_tolerance=0.0002
    mesh_args=(--step 0.05 --residual-tolerance "$vortex_residual_tolerance" --seed-grid 3 --max-triangles 100000)
    # (0.60,0), (0.90,0), (0.75,0.35), (0.75,-0.35), converted using
    # 2*(value-min)/(max-min)-1 independently on the u and v training binaries.
    edge_a=0.04409611104351585,-0.0001770810478360918
    edge_b=0.3659839633726256,-0.0001770810478360918
    wake_up=0.20504003720807096,0.39000094753311165
    wake_down=0.20504003720807096,-0.3903551096287836
}
