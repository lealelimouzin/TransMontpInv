######################################################################
# Leave the following line as it is
[ALL]
######################################################################

# DEFINE A NUMBER OF DIRECTORIES

# - directory containing files mksamples.in, procsamples.in,
#   forward_file, electrode_file, data_file
inversion_directory=~/erica/TDTOMO/ForTERT/bin_blocks

# - directory containing outputs from procsamples
file_directory=~/erica/TDTOMO/ForTERT/bin_blocks

# - directory containing forward modelling and interpolation codes
forward_directory=~/erica/TDTOMO/ForTERT/bin_blocks
gmtplot_directory=${forward_directory}/gmtplot

# - directory for output plots
plot_directory=~/erica/TDTOMO/ForTERT/examples/blocks/results/plots

######################################################################

# DEFINE A NUMBER OF PROGRAMS

# - program for map interpolation
map_interpolation=interpolate_2D

# - program for forward modelling
forward_modelling=forward_2_5d

# - program for calculation of the posterior
get_posterior=value_3D_post

######################################################################

# WHAT TO PLOT? ("yes" will plot)

# - average, median, maximum, harmonic mean and root-mean-square
plot_maps_1=yes

# - standard deviation, entropy, node density, skewness and kurtosis
plot_maps_2=yes

# - initial and final Voronoi models (one per chain)
plot_voronoi=yes

# - residuals
plot_residuals=yes

# - posterior on number of cells
plot_post_ncell=yes

# - posterior on noise parameters
plot_post_noise=yes

# - posterior on resistivities
plot_post_value=yes

# - number of cells and misfit vs. iteration
plot_chain_evol=yes

######################################################################

# GENERAL PARAMETERS FOR GMT SETUP

# - prefix for GMT commands
gmt_prefix=gmt

# - paper size (in GMT format for PS_MEDIA)
paper_size=Custom_30cx30c

# - thickness of the plot frame
frame_thickness=thin

# - position (in cm) of bottom left corner of maps/plots
x_corner=4
y_corner=3.4

# - format for output plot files (in GMT format)
plot_format="-Tg -E200"

######################################################################

# COLOR PALETTE PARAMETERS

# - location of color palette (.cpt) files
cpt_directory=~/erica/TDTOMO/ForTSWT/source/cptfiles

# - color palette, options, ticks (in GMT format) and colorbar label
#   for average, median, maximum, harmonic mean and root-mean-square
map_cpt=GMT_rainbow.cpt
map_cpt_options="-I -Z" # see -B flag in GMT commands
map_ticks=a1f0.5 # see -B flag in GMT commands
map_label="log@-10@-(@~r@~)"

# - color palette, options and ticks (in GMT format) for standard 
#   deviation and entropy
std_cpt=jet.cpt
std_cpt_options="-Z" # see -B flag in GMT commands
stdev_ticks=a0.5f0.1 # see -B flag in GMT commands
entropy_ticks=a1f0.5 # see -B flag in GMT commands

# - color palette and ticks (in GMT format) for node density
nod_cpt=GMT_hot.cpt
nod_cpt_options="-Z -I" # see -B flag in GMT commands
nod_ticks=a0.02f0.01 # see -B flag in GMT commands

# - color palette, options, ticks (in GMT format) and maximum value
#   of colorbar for skewness
skew_cpt=dkpurplegreen.cpt
skew_cpt_options="-Z"
max_abs_skew=10
skew_ticks=a5f1 # see -B flag in GMT commands

# - color palette, options, ticks (in GMT format) and maximum value
#   of colorbar for kurtosis
kurt_cpt=dkpurplegreen.cpt
kurt_cpt_options="-I"
kurt_ticks=a0.2+0.1 # see -B flag in GMT commands
kurt_sign_only=yes
zero_lim_kurt=0.1

######################################################################

# SETUP FOR MAPS

# - width and height of maps
map_width=6
map_height=0

# - x and y ticks (see -B flag in GMT commands)
map_x_ticks=a10f5
map_y_ticks=a10f5

# - x and y labels for maps (see -B flag in GMT commands)
map_x_label="x (m)"
map_y_label="z (m)"

# - position and height (in cm) of the color bar
x_cpt=7 #`echo ${x_corner} ${map_width} | awk '{print $1+$2/2}'`
y_cpt=2.6
cpt_height=0.2h

# - interpolation factor for maps and voronoi models
dicing_map=5
dicing_voro=10

# - plot electrodes? "yes" will plot them
plot_electrodes=yes
electrode_symbol="-Sy0.10 -Ggrey40 -W1.2,grey40 -N -D0/0.05" # symbol in GMT format

######################################################################

# SETUP FOR XY PLOTS

# - width and height of plots
xy_width=6
xy_height=4

# - x and y ticks for posterior on number of cells (in GMT format)
#   (if set to "", ticks will be assigned automatically)
ncell_x_ticks=a20f10
ncell_y_ticks=a0.02f0.01

# - x and y ticks for posterior on noise parameters (in GMT format)
#   (if set to "", ticks will be assigned automatically)
noise_x_ticks=a2f1
noise_y_ticks=a0.2f0.1

# - x and y ticks for actual residuals (in GMT format)
#   (if set to "", ticks will be assigned automatically)
resid_abs_x_ticks=a50f10
resid_abs_y_ticks=a1

# - x and y ticks for percentage residuals (in GMT format)
#   (if set to "", ticks will be assigned automatically)
resid_perc_x_ticks=a100f50
resid_perc_exp=4 # power of 10 for y axis of histogram
resid_perc_y_ticks=a2f0.5

######################################################################

# POINTS FOR POSTERIOR ON RESISTIVITY

# Choose points for which resistivity posteriors should be calculated
# and y ticks for plot in GMT format (if set to "", ticks will be 
# assigned automatically)
# (points assumed to be on the xz plane where y=0)

number_of_points=10

p1_x=-7.5
p1_y=0
p1_y_ticks=a4f2

p2_x=-7.5
p2_y=5
p2_y_ticks=a2f1

p3_x=-7.5
p3_y=10
p3_y_ticks=a1f0.5

p4_x=-7.5
p4_y=15
p4_y_ticks=a1f0.5

p5_x=-7.5
p5_y=20
p5_y_ticks=a1f0.5

p6_x=0
p6_y=0
p6_y_ticks=a2f1

p7_x=0
p7_y=5
p7_y_ticks=a2f1

p8_x=0
p8_y=10
p8_y_ticks=a2f1

p9_x=0
p9_y=15
p9_y_ticks=a1f0.5

p10_x=0
p10_y=20
p10_y_ticks=a1f0.5

######################################################################

# PARAMETERS FOR CHAIN EVOLUTION PLOTS

# How many independent Markov chains? (i.e. how many cores was the
# inversion run on in parallel?)
number_of_chains=6

# - x and y ticks for number of cells vs. iteration (in GMT format)
#   (if set to "", ticks will be assigned automatically)
chain_ncell_exp=3 # power of 10 for x axis
chain_ncell_x_ticks=a2f1
chain_ncell_y_ticks=a50f10

# - x and y ticks for number of cells vs. log10(misfit) (in GMT format)
#   (if set to "", ticks will be assigned automatically)
chain_misfits_exp=3 # power of 10 for x axis
chain_misfits_x_ticks=a2f1
chain_misfits_y_ticks=a1f0.5

######################################################################

# PARAMETERS FOR SYNTHETIC MODEL PLOTS

# Is this a synthetic test? ("yes" means it is)
synthetic_test=yes

# Do you want to plot the synthetic model? ("yes" means you do)
plot_synthetic=yes

# If it is a synthetic test, what is the real model file?
synthetic_directory=~/erica/TDTOMO/ForTERT/examples/blocks
synthetic_file=model3.vtx
synthetic_apparent_file=model3d_app.vtx

######################################################################

# PARAMETERS FOR LINEARISED INVERSION PLOTS

# Plot inversion results from R2? ("yes" means they are plotted)
plot_lin=yes

# What are the linearised inversion files?
lin_directory=~/erica/TDTOMO/ForTERT/examples/blocks/results/r2
lin_res_file=res.vtx
lin_unc_file=sen.vtx

# Label and ticks for colorbar of uncertainty map
lin_unc_label="log@-10@-(sensitivity)"
lin_unc_ticks=a1f0.5


######################################################################

# GET INFORMATION FROM mksamples.in AND procsamples.in

# Define path to mksamples.in and procsamples.in
mksamples_in=${inversion_directory}/mksamples.in
procsamples_in=${inversion_directory}/procsamples.in

# Input files for inversion (from mksamples.in)
electrode_file=`awk '{ if(NR==4) print $1}' ${mksamples_in}`
data_file=`awk '{ if(NR==5) print $1}' ${mksamples_in}`
modelling_method=`awk '{ if(NR==51) print $1}' ${mksamples_in}`
forward_file=`awk '{ if(NR==52) print $1}' ${mksamples_in}`

# Output files from inversion (from mksamples.in)
ncell_file=`awk '{ if(NR==12) print $1}' ${mksamples_in}`
misfits_file=`awk '{ if(NR==14) print $1}' ${mksamples_in}`

# Output files from inversion and processing (from procsamples.in)
average_file=`awk '{ if(NR==14) print $1}' ${procsamples_in}`
stdev_file=`awk '{ if(NR==15) print $1}' ${procsamples_in}`
median_file=`awk '{ if(NR==16) print $1}' ${procsamples_in}`
maximum_file=`awk '{ if(NR==17) print $1}' ${procsamples_in}`
harmean_file=`awk '{ if(NR==18) print $1}' ${procsamples_in}`
rms_file=`awk '{ if(NR==19) print $1}' ${procsamples_in}`
entropy_file=`awk '{ if(NR==20) print $1}' ${procsamples_in}`
skewness_file=`awk '{ if(NR==21) print $1}' ${procsamples_in}`
kurtosis_file=`awk '{ if(NR==22) print $1}' ${procsamples_in}`
node_file=`awk '{ if(NR==28) print $1}' ${procsamples_in}`
first_file=`awk '{ if(NR==29) print $1}' ${procsamples_in}`
last_file=`awk '{ if(NR==30) print $1}' ${procsamples_in}`
ncellsi_file=`awk '{ if(NR==23) print $1}' ${procsamples_in}`
misfiti_file=`awk '{ if(NR==24) print $1}' ${procsamples_in}`
postncells_file=`awk '{ if(NR==25) print $1}' ${procsamples_in}`
postnoise_file=`awk '{ if(NR==26) print $1}' ${procsamples_in}`
postvalue_file=`awk '{ if(NR==27) print $1}' ${procsamples_in}`
residuals_histo_file=`awk '{ if(NR==31) print $1}' ${procsamples_in}`

# Prior (from mksamples.in)
x_min=`awk '{ if(NR==24) print $1}' ${mksamples_in}`
x_max=`awk '{ if(NR==24) print $2}' ${mksamples_in}`
y_min=`awk '{ if(NR==25) print $1}' ${mksamples_in}`
y_max=`awk '{ if(NR==25) print $2}' ${mksamples_in}`
z_min=`awk '{ if(NR==26) print $1}' ${mksamples_in}`
z_max=`awk '{ if(NR==26) print $2}' ${mksamples_in}`
ncell_min=`awk '{ if(NR==22) print $1}' ${mksamples_in}`
ncell_max=`awk '{ if(NR==22) print $2}' ${mksamples_in}`
value_min=`awk '{ if(NR==23) print $1}' ${mksamples_in}`
value_max=`awk '{ if(NR==23) print $2}' ${mksamples_in}`
number_of_datasets=`awk '{ if(NR==6) print $1}' ${mksamples_in}`
noise_type=`awk '{ if(NR==27) print $1}' ${mksamples_in}`
a_min=`awk '{ if(NR==28) print $1}' ${mksamples_in}`
a_max=`awk '{ if(NR==28) print $2}' ${mksamples_in}`
b_min=`awk '{ if(NR==29) print $1}' ${mksamples_in}`
b_max=`awk '{ if(NR==29) print $2}' ${mksamples_in}`
lambda_min=`awk '{ if(NR==30) print $1}' ${mksamples_in}`
lambda_max=`awk '{ if(NR==30) print $2}' ${mksamples_in}`
sigma_min=`awk '{ if(NR==31) print $1}' ${mksamples_in}`
sigma_max=`awk '{ if(NR==31) print $2}' ${mksamples_in}`

# Map pixelization parameters
x_number_of_pixels=`awk '{ if(NR==36) print $1}' ${procsamples_in}`
y_number_of_pixels=`awk '{ if(NR==36) print $2}' ${procsamples_in}`
z_number_of_pixels=`awk '{ if(NR==36) print $3}' ${procsamples_in}`
