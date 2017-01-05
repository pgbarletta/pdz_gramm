# --------------------------------------------
# Copy ligand bound files to post qha analysis
# --------------------------------------------

# Copy ligand bound qha output
cp ../qha/lb_prod_vecs_ca.dat .
# Copy ligand bound qha covar matrix
cp ../qha/lb_prod_covar_ca.dat .
# Copy ligand bound traj to qha post analysis
cp ../qha/lb_prod_r_ca.mdcrd .
# Copy ligand bound topology and change name to qha post analysis
cp /scratch/lfs/grosso/gip/parmtop_files/lf_strip_tails.top lb_prod.top
