# --------------------------------------------
# Copy ligand free files to post qha analysis
# --------------------------------------------

# Copy ligand free qha output
cp ../qha/lf_prod_vecs_ca.dat .
# Copy ligand free qha covar matrix
cp ../qha/lf_prod_covar_ca.dat .
# Copy ligand free traj to qha post analysis
cp ../qha/lf_prod_r_ca.mdcrd .
# Copy ligand free topology and change name to qha post analysis
cp /scratch/lfs/grosso/gip/parmtop_files/lf_strip_tails.top lf_prod.top
