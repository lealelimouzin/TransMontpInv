#!/bin/bash

# ********************************************************************* #
# Script to plot outputs from procsamples.f90 using GMT
# Copyright (C) 2017  Erica Galetti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ********************************************************************* #


# Get parameters script
if [[ ${1} == "" ]] ; then
        echo
        echo "   Argument of ./plot_TERT.sh must be"
        echo "   the file of plotting parameters!"
        echo "   TERMINATING PROGRAM!!"
        echo
        exit
fi
sed -i 's/ = /=/g' ${1} # ensure there are no spaces around "=" signs
. <(grep -v "^\[" ${1}) # execute file of parameters (removing section headings)

# Output date and time
echo
echo "*************************************************"
echo "Plotting started:"
date "+%d/%m/%Y %H:%M:%S"
echo "*************************************************"
echo

# Electrode coordinats are given as [x y].
xy_toggle=""

# If map_height=0 (i.e., map height is computed from the map width),
# calculate the height of the map and set it positive downward
if [[ ${map_height} == "0" ]] ; then
        map_height=`echo ${x_min} ${x_max} ${z_min} ${z_max} ${map_width} | awk '{print $5*($4-$3)/($2-$1)}'`
fi

# Set projections:
# - type of projection for maps (in GMT format)
projection_map=-JX${map_width}/-${map_height}
# - type of projection for xy plots (i.e., histograms and chain 
#   evolution) (in GMT format)
projection_xy=-JX${xy_width}/${xy_height}

# Create directory for outputs
mkdir -p ${plot_directory}

# Move to gmt-plotting directory and set up GMT
cd ${gmtplot_directory}
${gmt_prefix} gmtset FORMAT_GEO_MAP ddd \
                     PROJ_LENGTH_UNIT cm \
                     MAP_FRAME_TYPE plain \
                     MAP_FRAME_PEN ${frame_thickness},black \
                     MAP_FRAME_AXES WeSn \
                     FONT_LABEL 6 \
                     FONT_ANNOT_PRIMARY 6 \
                     MAP_LABEL_OFFSET 0.15c \
                     MAP_TICK_LENGTH -0.1c \
                     PS_PAGE_ORIENTATION portrait \
                     PS_MEDIA ${paper_size}
psfile=plot.ps # name of ps file
pngfile=plot.png # name of png file
color_increment=0.001

# Move to forward-modelling directory and set up forward modelling
cd ${forward_directory}
cp ${inversion_directory}/${electrode_file} electrodes.dat
cp ${inversion_directory}/${data_file} data.dat
cp ${inversion_directory}/${forward_file} forward.in
# Set up forward_2_5d.in
awk '{ if(NR==4) print "electrodes.dat                :: File of electrode coordinates";
  else if(NR==5) print "data.dat                      :: File of resistance measurements";
  else if(NR==6) print "1                             :: Model file type (0=smooth, 1=voronoi)";
  else if(NR==7) print "gridc.vtx                     :: Resistivity model file";
  else if(NR==8) print "1                             :: Resistivities in model file are log10 (1=yes, 0=no)";
  else if(NR==9) print "forward.in                    :: Forward modelling file";
  else if(NR==13) print "resistances.out               :: File of resistances";
  else if(NR==14) print "apparentrho.out               :: File of apparent resistivities";
  else print $0 }
  ' ${forward_modelling}.in > ${forward_modelling}1.in
mv ${forward_modelling}1.in ${forward_modelling}.in

# Move to gmt-plotting directory and set up interpolate_2D.in, with interpolation parameters
cd ${gmtplot_directory}
rm -f *.cpt *.ps
# Set up interpolate_2D.in
awk -v npx=`echo ${x_number_of_pixels}` -v npy=`echo ${y_number_of_pixels}` -v npz=`echo ${z_number_of_pixels}` \
    -v dic=`echo ${dicing_map}` \
    -v xmin=`echo ${x_min}` -v xmax=`echo ${x_max}` \
    -v ymin=`echo ${y_min}` -v ymax=`echo ${y_max}` \
    -v zmin=`echo ${z_min}` -v zmax=`echo ${z_max}` '
        { if(NR==4) print "2                              :: Model file format: Voronoi (=1), vtx 2D (=2), vtx 3D (=3)";
     else if(NR==5) print "../gridc.vtx                      :: Model file";
     else if(NR==6) print "../electrodes.dat              :: Electrodes file";
     else if(NR==10) print dic "                             :: Grid dicing";
     else if(NR==11) print "map.z                          :: GMT output file for slice";
     else if(NR==12) print "bound.gmt                      :: GMT plotting bounds for slice";
     else if(NR==13) print "1                              :: Plot electrodes?";
     else if(NR==14) print "electrodes.dat                 :: Electrodes file compatible with GMT";
     else if(NR==18) print "1                              :: Plot xy plane (=0), xz plane (=1), yz plane (=2)";
     else if(NR==19) print "0                              :: Value of z if xy plane, y if xz plane, x if yz plane";
     else if(NR==23) print xmin "   " xmax "                       :: Minimum and maximum x";
     else if(NR==24) print ymin "   " ymax "                       :: Minimum and maximum y";
     else if(NR==25) print zmin "   " zmax "                       :: Minimum and maximum z";
     else if(NR==26) print npx "   " npy "   " npz "                   :: Number of pixels in [x y z]";
     else print $0 }
     ' ${map_interpolation}.in > ${map_interpolation}1.in
mv ${map_interpolation}1.in ${map_interpolation}.in

# Print message to screen
echo "-------------------- PLOTTING -------------------"

# Plot average, median, maximum, harmonic mean and root-mean-square
if [[ ${plot_maps_1} == "yes" ]] ; then
        for map_file in ${average_file} ${median_file} ${maximum_file} ${harmean_file} ${rms_file} ${synthetic_file} ${synthetic_apparent_file} ${lin_res_file} ; do
                
                if [[ ${synthetic_test} != "yes" ]] || [[ ${plot_synthetic} != "yes" ]] ; then
                        if [[ ${map_file} == "${synthetic_file}" ]] || [[ ${map_file} == "${synthetic_apparent_file}" ]] ; then
                                continue
                        fi
                fi
                if [[ ${plot_lin} != "yes" ]] ; then
                        if [[ ${map_file} == "${lin_res_file}" ]] ; then
                                continue
                        fi
                fi
                
                # Set up names of image files
                map_image=`echo "${map_file}" | awk 'BEGIN {FS="."} ; {print $1}'`
                
                # Move to gmt-plotting directory to plot
                cd ${gmtplot_directory}
                
                # Copy grid file and set up plotting parameters
                if [[ ${map_file} == "${synthetic_file}" ]] ; then
                        
                        cp ${synthetic_directory}/${synthetic_file} ../gridc.vtx
                        
                        awk '{if(NR==10) print "1                             :: Grid dicing";
                              else print $0 }
                            ' ${map_interpolation}.in > ${map_interpolation}1.in
                            mv ${map_interpolation}1.in ${map_interpolation}.in
                            
                elif [[ ${map_file} == "${synthetic_apparent_file}" ]] ; then
                        
                        cp ${synthetic_directory}/${synthetic_apparent_file} ../gridc.vtx
                        
                        awk '{if(NR==10) print "1                             :: Grid dicing";
                              else print $0 }
                            ' ${map_interpolation}.in > ${map_interpolation}1.in
                            mv ${map_interpolation}1.in ${map_interpolation}.in
                            
                else
                        
                        if [[ ${map_file} == "${lin_res_file}" ]] ; then
                                cp ${lin_directory}/${lin_res_file} ../gridc.vtx
                        else
                                cp ${file_directory}/${map_file} ../gridc.vtx
                        fi
                        awk -v dic=`echo ${dicing_map}` \
                                '{if(NR==10) print dic "                             :: Grid dicing";
                                  else print $0 }
                                ' ${map_interpolation}.in > ${map_interpolation}1.in
                                mv ${map_interpolation}1.in ${map_interpolation}.in
                                
                fi
                
                # Set up color palette label
                if [[ ${map_file} == "${average_file}" ]] ; then
                        label="arithmetic mean of ${map_label}"
                elif [[ ${map_file} == "${median_file}" ]] ; then
                        label="median of ${map_label}"
                elif [[ ${map_file} == "${maximum_file}" ]] ; then
                        label="mode of ${map_label}"
                elif [[ ${map_file} == "${harmean_file}" ]] ; then
                        label="harmonic mean of ${map_label}"
                elif [[ ${map_file} == "${rms_file}" ]] ; then
                        label="root-mean-square of ${map_label}"
                elif [[ ${map_file} == "${synthetic_apparent_file}" ]] ; then
                        label="log@-10@-(@~r@~@-app@-)"
                else
                        label=${map_label}
                fi
                
                # Get interpolated map
                ./${map_interpolation}
                
                # Get model boundaries
                # - boundary coords from bound.gmt
                left=`awk 'NR==1 {printf "%.2f",$1}' bound.gmt`
                right=`awk 'NR==2 {printf "%.2f",$1}' bound.gmt`
                bottom=`awk 'NR==3 {printf "%.2f",$1}' bound.gmt`
                top=`awk 'NR==4 {printf "%.2f",$1}' bound.gmt`
                # - gridding units
                x_unit=`awk 'NR==5 {print $1}' bound.gmt`
                y_unit=`awk 'NR==6 {print $1}' bound.gmt`
                # Projection boundaries
                if [[ ${map_file} == "${synthetic_file}" ]] || [[ ${map_file} == "${synthetic_apparent_file}" ]] || [[ ${map_file} == "${lin_res_file}" ]] ; then
                        boundaries_ext=-R${left}/${right}/${bottom}/${top}
                else
                        boundaries=-R${left}/${right}/${bottom}/${top}
                fi
                
                # Set color palette boundaries
                up=`echo ${value_max} + ${color_increment} | bc -l`
                down=`echo ${value_min} - ${color_increment} | bc -l`
                # Create colour palette
                cp ${cpt_directory}/${map_cpt} .
                cpt=${map_cpt}
                ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${map_cpt_options} > ${gmtplot_directory}/palette.cpt
                
                # Plot the map
                echo "-------------------------------------------------"
                echo "Plotting ${map_file}..."
                
                if [[ ${map_file} == "${synthetic_file}" ]] || [[ ${map_file} == "${synthetic_apparent_file}" ]] || [[ ${map_file} == "${lin_res_file}" ]] ; then
                        ${gmt_prefix} xyz2grd map.z -Gmap.grd -I${x_unit}/${y_unit} -ZLB ${boundaries_ext}
                else
                        ${gmt_prefix} xyz2grd map.z -Gmap.grd -I${x_unit}/${y_unit} -ZLB ${boundaries}
                fi
                ${gmt_prefix} grdimage map.grd ${boundaries} ${projection_map} -Cpalette.cpt -Q -K -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                ${gmt_prefix} psbasemap ${boundaries} ${projection_map} -B${map_x_ticks}:"${map_x_label}":/${map_y_ticks}:"${map_y_label}": -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                if [[ ${plot_electrodes} == "yes" ]] ; then
                        ${gmt_prefix} psxy electrodes.dat ${boundaries} ${projection_map} ${xy_toggle} ${electrode_symbol} -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                fi
                # Plot color scale
                ${gmt_prefix} psscale -Cpalette.cpt -B${map_ticks}:"${label}": -D${x_cpt}/${y_cpt}/${map_width}/${cpt_height} -O -P >> ${psfile}
                # Convert ps to png and move file
                ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${map_image}
                
                echo "...done!"
                
        done
fi

# Plot standard deviation, entropy, node density, skewness and kurtosis
if [[ ${plot_maps_2} == "yes" ]] ; then
        for map_file in ${stdev_file} ${entropy_file} ${node_file} ${skewness_file} ${kurtosis_file} ${lin_unc_file} ; do
        #for map_file in ${skewness_file} ${kurtosis_file} ; do
        #for map_file in ${stdev_file} ${lin_unc_file} ; do
                
                if [[ ${plot_lin} != "yes" ]] ; then
                        if [[ ${map_file} == "${lin_unc_file}" ]] ; then
                                continue
                        fi
                fi
                # Set up names of image files
                map_image=`echo "${map_file}" | awk 'BEGIN {FS="."} ; {print $1}'`
                
                # Move to gmt-plotting directory to plot
                cd ${gmtplot_directory}
                        
                # Copy grid file
                if [[ ${map_file} == "${lin_unc_file}" ]] ; then
                        cp ${lin_directory}/${lin_unc_file} ../gridc.vtx
                else
                        cp ${file_directory}/${map_file} ../gridc.vtx
                fi
                
                # Get interpolated map
                ./${map_interpolation}
                
                # Get model boundaries
                # - boundary coords from bound.gmt
                left=`awk 'NR==1 {printf "%.2f",$1}' bound.gmt`
                right=`awk 'NR==2 {printf "%.2f",$1}' bound.gmt`
                bottom=`awk 'NR==3 {printf "%.2f",$1}' bound.gmt`
                top=`awk 'NR==4 {printf "%.2f",$1}' bound.gmt`
                # - gridding units
                x_unit=`awk 'NR==5 {print $1}' bound.gmt`
                y_unit=`awk 'NR==6 {print $1}' bound.gmt`
                # Projection boundaries
                if [[ ${map_file} == "${lin_unc_file}" ]] ; then
                        boundaries_ext=-R${left}/${right}/${bottom}/${top}
                else
                        boundaries=-R${left}/${right}/${bottom}/${top}
                fi
                
                # Set up color bar
                grep -v "NaN" map.z > tmp.z
                if [[ ${map_file} == "${stdev_file}" ]] || [[ ${map_file} == "${entropy_file}" ]] ; then
                        # Create colour palette
                        cp ${cpt_directory}/${std_cpt} .
                        cpt=${std_cpt}
                        if [[ ${map_file} == "${stdev_file}" ]] ; then
                                up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max+'${color_increment}'}' tmp.z`                
                                down=0
                                ticks=${stdev_ticks}
                                label="standard deviation of ${map_label}"
                                ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${std_cpt_options} > ${gmtplot_directory}/stdev_palette.cpt
                        elif [[ ${map_file} == "${entropy_file}" ]] ; then
                                up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max+'${color_increment}'}' tmp.z`                
                                down=0 #`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min-'${color_increment}'}' tmp.z`
                                ticks=${entropy_ticks}
                                label="entropy of ${map_label}"
                        fi
                        ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${std_cpt_options} > ${gmtplot_directory}/palette.cpt
                elif [[ ${map_file} == "${node_file}" ]] ; then
                        # Create colour palette
                        cp ${cpt_directory}/${nod_cpt} .
                        cpt=${nod_cpt}
                        up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max+'${color_increment}'}' tmp.z`                
                        down=0 #`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min-'${color_increment}'}' tmp.z`
                        ticks=${nod_ticks}
                        label="node density"
                        ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${nod_cpt_options} > ${gmtplot_directory}/palette.cpt
                elif [[ ${map_file} == "${skewness_file}" ]] || [[ ${map_file} == "${kurtosis_file}" ]] ; then
                        # Create colour palette
                        if [[ ${map_file} == "${skewness_file}" ]] ; then
                                cp ${cpt_directory}/${skew_cpt} .
                                cpt=${skew_cpt}
                                if [[ ${max_abs_skew} != "" ]] ; then
                                        awk -v sm=`echo ${max_abs_skew}` \
                                           '{if ($1<-sm) print -sm;
                                              else if ($1>sm) print sm;
                                              else print $0}' map.z > tmp.z
                                        mv tmp.z map.z
                                fi
                                grep -v "NaN" map.z > tmp.z
                                up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max+'${color_increment}'}' tmp.z`                
                                down=`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min-'${color_increment}'}' tmp.z`
                                echo ${up} | awk '{ if($1>0) print $1; else print -1*$1}' > tmp.txt
                                echo ${down} | awk '{ if($1>0) print $1; else print -1*$1}' >> tmp.txt
                                up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max}' tmp.txt`
                                down=-${up}
                                rm -f tmp.txt
                                ticks=${skew_ticks}
                                label="skewness of ${map_label}"
                                ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${skew_cpt_options} > ${gmtplot_directory}/palette.cpt
                        elif [[ ${map_file} == "${kurtosis_file}" ]] ; then
                                cp ${cpt_directory}/${kurt_cpt} .
                                cpt=${kurt_cpt}
                                if [[ ${kurt_sign_only} == "yes" ]] ; then
                                        awk -v zero=`echo ${zero_lim_kurt}` \
                                           '{if ($1<-zero) print -zero*2;
                                              else if ($1>zero) print zero*2;
                                              else if ($1>=-zero && $1<=zero) print 0;
                                              else print $0}' map.z > tmp.z
                                        mv tmp.z map.z
                                else
                                        awk '{if ($1<=0) print NaN;
                                              else print log($1+3)/log(10)}' map.z > tmp.z
                                        mv tmp.z map.z
                                fi
                                grep -v "NaN" map.z > tmp.z
                                if [[ ${kurt_sign_only} == "yes" ]] ; then
                                        up=`echo ${zero_lim_kurt} | awk '{print $1*3}'`
                                        down=-${up}
                                        inc=`echo ${zero_lim_kurt} | awk '{print $1*2}'`
                                        ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${inc} ${kurt_cpt_options} > ${gmtplot_directory}/palette.cpt
                                        label="excess kurtosis of ${map_label}"
                                else
                                        up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max+'${color_increment}'}' tmp.z`                
                                        down=0 #`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min-'${color_increment}'}' tmp.z`
                                        #echo ${up} | awk '{ if($1>0) print $1; else print -1*$1}' > tmp.txt
                                        #echo ${down} | awk '{ if($1>0) print $1; else print -1*$1}' >> tmp.txt
                                        #up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max}' tmp.txt`
                                        #down=-${up}
                                        ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${kurt_cpt_options} > ${gmtplot_directory}/palette.cpt
                                        label="kurtosis"
                                fi
                                rm -f tmp.txt
                                ticks=${kurt_ticks}
                        fi
                elif [[ ${map_file} == "${lin_unc_file}" ]] ; then
                        # Create colour palette
                        cp ${cpt_directory}/${std_cpt} .
                        cpt=${std_cpt}
                        #up=0
                        up=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max+'${color_increment}'}' tmp.z`
                        down=`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min-'${color_increment}'}' tmp.z`
                        ticks=${lin_unc_ticks}
                        label=${lin_unc_label}
                        ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} ${std_cpt_options} -I > ${gmtplot_directory}/palette.cpt
                fi
                
                # Plot the map
                echo "-------------------------------------------------"
                echo "Plotting ${map_file}..."
                
                if [[ ${map_file} == "${lin_unc_file}" ]] ; then
                        ${gmt_prefix} xyz2grd map.z -Gmap.grd -I${x_unit}/${y_unit} -ZLB ${boundaries_ext}
                else
                        ${gmt_prefix} xyz2grd map.z -Gmap.grd -I${x_unit}/${y_unit} -ZLB ${boundaries}
                fi
                ${gmt_prefix} grdimage map.grd ${boundaries} ${projection_map} -Cpalette.cpt -Q -K -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                ${gmt_prefix} psbasemap ${boundaries} ${projection_map} -B${map_x_ticks}:"${map_x_label}":/${map_y_ticks}:"${map_y_label}": -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                if [[ ${plot_electrodes} == "yes" ]] ; then
                        ${gmt_prefix} psxy electrodes.dat ${boundaries} ${projection_map} ${xy_toggle} ${electrode_symbol} -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                fi
                # Plot color scale
                ${gmt_prefix} psscale -Cpalette.cpt -B${ticks}:"${label}": -D${x_cpt}/${y_cpt}/${map_width}/${cpt_height} -O -P >> ${psfile}
                # Convert ps to png and move file
                ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${map_image}
                
                echo "...done!"
                
        done
fi

# Plot first and last model on each chain
if [[ ${plot_voronoi} == "yes" ]] ; then
        for chain_number in $(seq 1 1 ${number_of_chains}) ; do
                
                # Move to gmt-plotting directory to plot
                cd ${gmtplot_directory}
                
                # Get interpolated map for first model
                awk -v dic=`echo ${dicing_voro}` '
                   { if(NR==4) print "1                              :: Model file format: Voronoi (=1), vtx 2D (=2), vtx 3D (=3)";
                else if(NR==5) print "../gridc.vor                      :: Model file";
                else if(NR==10) print dic "                             :: Grid dicing";
                else print $0 }
                ' ${map_interpolation}.in > ${map_interpolation}1.in
                mv ${map_interpolation}1.in ${map_interpolation}.in
                
                echo "-------------------------------------------------"
                echo "Plotting initial and final model, chain ${chain_number}..."
                map_image=first_last_${chain_number}
                
                # Initial model
                cp ${file_directory}/${first_file}${chain_number} ../gridc.vor
                ./${map_interpolation}
                mv map.z mapf.z
                
                # Final model
                cp ${file_directory}/${last_file}${chain_number} ../gridc.vor
                ./${map_interpolation}
                mv map.z mapl.z
                                
                # Get model boundaries
                # - boundary coords from bound.gmt
                left=`awk 'NR==1 {printf "%.2f",$1}' bound.gmt`
                right=`awk 'NR==2 {printf "%.2f",$1}' bound.gmt`
                bottom=`awk 'NR==3 {printf "%.2f",$1}' bound.gmt`
                top=`awk 'NR==4 {printf "%.2f",$1}' bound.gmt`
                # - gridding units
                x_unit=`awk 'NR==5 {print $1}' bound.gmt`
                y_unit=`awk 'NR==6 {print $1}' bound.gmt`
                # Projection boundaries
                boundaries=-R${left}/${right}/${bottom}/${top}
                
                # Set color palette boundaries
                up=`echo ${value_max} + ${color_increment} | bc -l`
                down=`echo ${value_min} - ${color_increment} | bc -l`
                # Create colour palette
                cp ${cpt_directory}/${map_cpt} .
                cpt=${map_cpt}
                ${gmt_prefix} makecpt -C${cpt} -T${down}/${up}/${color_increment} -I -Z > ${gmtplot_directory}/palette.cpt
                
                # Plot first model
                spacing=0.5
                y_corner2=`echo ${y_corner} ${map_height} ${spacing} | awk '{print $1+$2+$3}'`
                ${gmt_prefix} gmtset MAP_FRAME_AXES Wesn
                ${gmt_prefix} xyz2grd mapf.z -Gmap.grd -I${x_unit}/${y_unit} -ZLB ${boundaries}
                ${gmt_prefix} grdimage map.grd ${boundaries} ${projection_map} -Cpalette.cpt -Q -K -P -Xa${x_corner} -Ya${y_corner2} > ${psfile}
                ${gmt_prefix} psbasemap ${boundaries} ${projection_map} -B${map_x_ticks}:"${map_x_label}":/${map_y_ticks}:"${map_y_label}": -O -K -P -Xa${x_corner} -Ya${y_corner2} >> ${psfile}
                if [[ ${plot_electrodes} == "yes" ]] ; then
                        ${gmt_prefix} psxy electrodes.dat ${boundaries} ${projection_map} ${xy_toggle} ${electrode_symbol} -O -K -P -Xa${x_corner} -Ya${y_corner2} >> ${psfile}
                fi
                # Plot last model
                ${gmt_prefix} gmtset MAP_FRAME_AXES WeSn
                ${gmt_prefix} xyz2grd mapl.z -Gmap.grd -I${x_unit}/${y_unit} -ZLB ${boundaries}
                ${gmt_prefix} grdimage map.grd ${boundaries} ${projection_map} -Cpalette.cpt -Q -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                ${gmt_prefix} psbasemap ${boundaries} ${projection_map} -B${map_x_ticks}:"${map_x_label}":/${map_y_ticks}:"${map_y_label}": -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                if [[ ${plot_electrodes} == "yes" ]] ; then
                        ${gmt_prefix} psxy electrodes.dat ${boundaries} ${projection_map} ${xy_toggle} ${electrode_symbol} -O -K -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                fi
                # Plot color scale
                ${gmt_prefix} psscale -Cpalette.cpt -B${map_ticks}:"${map_label}": -D${x_cpt}/${y_cpt}/${map_width}/${cpt_height} -O -P >> ${psfile}
                # Convert ps to png and move file
                ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${map_image}
                
                echo "...done!"
                
        done
fi
rm -f map?.z

# Move to gmt-plotting directory to plot
cd ${gmtplot_directory}

${gmt_prefix} gmtset MAP_FRAME_AXES WeSn

# Plot residuals histogram
if [[ ${plot_residuals} == "yes" ]] ; then
        
        echo "-------------------------------------------------"
        echo "Plotting residuals..."
        
        awk '{print $1,log($2)/log(10),$3,$4/(10^'${resid_perc_exp}')}' ${file_directory}/${residuals_histo_file} > histo.dat
        # Absolute residuals
        histo_image=`echo "${residuals_histo_file}" | awk 'BEGIN {FS="."} ; {print $1"_abs"}'`
        xmin=`awk 'min=="" || $1 < min {min=$1} END {printf "%.4f",min}' histo.dat`
        xmax=`awk 'max=="" || $1 > max {max=$1} END {printf "%.4f",max}' histo.dat`
        ymax=`awk 'max=="" || $2 > max {max=$2} END {printf "%.4f",max}' histo.dat`
        boundaries=-R${xmin}/${xmax}/0/${ymax}
        if [[ ${resid_abs_x_ticks} == "" ]] ; then
                inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.2f",($2-$1)/4}'`
                inc_x2=`echo "${inc_x}" | awk '{printf "%.2f",($1)/2}'`
                resid_abs_x_ticks=a${inc_x}f${inc_x2}
        fi
        if [[ ${resid_abs_y_ticks} == "" ]] ; then
                inc_y=`echo "${ymax}" | awk '{printf "%d",($1)/4}'`
                inc_y2=`echo "${inc_y}" | awk '{printf "%.2f",($1)/2}'`
                resid_abs_y_ticks=a${inc_y}f${inc_y2}
        fi
        bin_width=`awk '{if (NR<=2) print $1}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
        awk '{print $1,$2}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                           -B${resid_abs_x_ticks}:"residual (@~W@~)":/${resid_abs_y_ticks}:"log@-10@-(number of readings)": \
                                                           -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
        # Convert ps to png and move file
        ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
        # Residuals as percentage of observed traveltimes
        histo_image=`echo "${residuals_histo_file}" | awk 'BEGIN {FS="."} ; {print $1"_perc"}'`
        xmin=`awk 'min=="" || $3 < min {min=$3} END {printf "%.4f",min}' histo.dat`
        xmax=`awk 'max=="" || $3 > max {max=$3} END {printf "%.4f",max}' histo.dat`
        ymax=`awk 'max=="" || $4 > max {max=$4} END {printf "%.4f",max}' histo.dat`
        boundaries=-R${xmin}/${xmax}/0/${ymax}
        if [[ ${resid_abs_x_ticks} == "" ]] ; then
                inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.2f",($2-$1)/4}'`
                inc_x2=`echo "${inc_x}" | awk '{printf "%.2f",($1)/2}'`
                resid_perc_x_ticks=a${inc_x}f${inc_x2}
        fi
        if [[ ${resid_abs_y_ticks} == "" ]] ; then
                inc_y=`echo "${ymax}" | awk '{printf "%d",($1)/4}'`
                inc_y2=`echo "${inc_y}" | awk '{printf "%.2f",($1)/2}'`
                resid_perc_y_ticks=a${inc_y}f${inc_y2}
        fi
        bin_width=`awk '{if (NR<=2) print $3}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
        awk '{print $3,$4}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                -B${resid_perc_x_ticks}:"residual (%)":/${resid_perc_y_ticks}:"number of readings (x10@+${resid_perc_exp}@+)": \
                                                                -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
        # Convert ps to png and move file
        ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
        
        echo "...done!"
        
fi

# Plot histogram of posterior on number of cells
if [[ ${plot_post_ncell} == "yes" ]] ; then
        
        echo "-------------------------------------------------"
        echo "Plotting posterior on number of cells..."
        
        cp ${file_directory}/${postncells_file} histo.dat
        histo_image=`echo "${postncells_file}" | awk 'BEGIN {FS="."} ; {print $1}'`
        xmin=${ncell_min}
        xmax=${ncell_max}
        ymax=`awk 'max=="" || $2 > max {max=$2} END {printf "%.6f",max}' histo.dat`
        boundaries=-R${xmin}/${xmax}/0/${ymax}
        if [[ ${ncell_x_ticks} == "" ]] ; then
                inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%d",($2-$1)/4}'`
                inc_x2=`echo "${inc_x}" | awk '{printf "%.2f",($1)/2}'`
                ncell_x_ticks=a${inc_x}f${inc_x2}
        fi
        if [[ ${ncell_y_ticks} == "" ]] ; then
                inc_y=`echo "${ymax}" | awk '{printf "%d",($1)/4}'`
                inc_y2=`echo "${inc_y}" | awk '{printf "%.2f",($1)/2}'`
                ncell_y_ticks=a${inc_y}f${inc_y2}
        fi
        awk '{print $1,$2}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                           -B${ncell_x_ticks}:"number of cells":/${ncell_y_ticks}:"p(n|d)": \
                                                           -Sb1u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
        # Convert ps to png and move file
        ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
        
        echo "...done!"
        
fi

# Plot histogram(s) of posterior on data noise
if [[ ${plot_post_noise} == "yes" ]] ; then
        
        cp ${file_directory}/${postnoise_file} histo.dat
        
        for dataset in $(seq 1 1 ${number_of_datasets})  ; do
                
                if [[ ${noise_type} == "1" || ${noise_type} == "2" ]] ; then
                        
                        echo "-------------------------------------------------"
                        echo "Plotting posterior on noise parameter a..."
                        
                        histo_image=`echo "${postnoise_file}" | awk 'BEGIN {FS="."} ; {print $1"_a_d'${dataset}'"}'`
                        xmin=${a_min}
                        xmax=${a_max}
                        column_x=1
                        column_y=`echo "${dataset}" | awk '{print $1+1}'`
                        ymax=`awk 'max=="" || $'${column_y}' > max {max=$'${column_y}'} END {printf "%.6f",max}' histo.dat`
                        boundaries=-R${xmin}/${xmax}/0/${ymax}
                        
                        if [[ ${noise_x_ticks} == "" ]] ; then
                                inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.4f",($2-$1)/4}'`
                                inc_x2=`echo "${inc_x}" | awk '{printf "%.4f",($1)/2}'`
                                noise_x_ticks=a${inc_x}f${inc_x2}
                        fi
                        if [[ ${noise_y_ticks} == "" ]] ; then
                                inc_y=`echo "${ymax}" | awk '{printf "%.4f",($1)/4}'`
                                inc_y2=`echo "${inc_y}" | awk '{printf "%.4f",($1)/2}'`
                                noise_y_ticks=a${inc_y}f${inc_y2}
                        fi
                        
                        bin_width=`awk '{if (NR<=2) print $'${column_x}'}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
                        awk '{print $'${column_x}',$'${column_y}'}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                                                   -B${noise_x_ticks}:"a":/${noise_y_ticks}:"p(a|d)": \
                                                                                                   -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                        # Convert ps to png and move file
                        ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
                        
                        echo "...done!"
                        
                        echo "-------------------------------------------------"
                        echo "Plotting posterior on noise parameter b..."
                        
                        histo_image=`echo "${postnoise_file}" | awk 'BEGIN {FS="."} ; {print $1"_b_d'${dataset}'"}'`
                        xmin=${b_min}
                        xmax=${b_max}
                        column_x=`echo "${number_of_datasets}" | awk '{print $1+2}'`
                        column_y=`echo "${dataset}" "${number_of_datasets}" | awk '{print $1+2+$2}'`
                        ymax=`awk 'max=="" || $'${column_y}' > max {max=$'${column_y}'} END {printf "%.6f",max}' histo.dat`
                        boundaries=-R${xmin}/${xmax}/0/${ymax}
                        
                        if [[ ${noise_x_ticks} == "" ]] ; then
                                inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.4f",($2-$1)/4}'`
                                inc_x2=`echo "${inc_x}" | awk '{printf "%.4f",($1)/2}'`
                                noise_x_ticks=a${inc_x}f${inc_x2}
                        fi
                        if [[ ${noise_y_ticks} == "" ]] ; then
                                inc_y=`echo "${ymax}" | awk '{printf "%.4f",($1)/4}'`
                                inc_y2=`echo "${inc_y}" | awk '{printf "%.4f",($1)/2}'`
                                noise_y_ticks=a${inc_y}f${inc_y2}
                        fi
                        
                        bin_width=`awk '{if (NR<=2) print $'${column_x}'}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
                        awk '{print $'${column_x}',$'${column_y}'}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                                                   -B${noise_x_ticks}:"b":/${noise_y_ticks}:"p(b|d)": \
                                                                                                   -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                        # Convert ps to png and move file
                        ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
                        
                        echo "...done!"
                        
                elif [[ ${noise_type} == "3" || ${noise_type} == "0" ]] ; then
                
                        if [[ ${noise_type} == "3" ]] ; then
                                
                                echo "-------------------------------------------------"
                                echo "Plotting posterior on noise parameter lambda..."
                                
                                xmin=${lambda_min}
                                xmax=${lambda_max}
                        elif [[ ${noise_type} == "0" ]] ; then
                                
                                echo "-------------------------------------------------"
                                echo "Plotting posterior on noise parameter sigma..."
                                
                                xmin=${sigma_min}
                                xmax=${sigma_max}
                        fi
                        
                        histo_image=`echo "${postnoise_file}" | awk 'BEGIN {FS="."} ; {print $1"_d'${dataset}'"}'`
                        column_x=1
                        column_y=`echo "${dataset}" | awk '{print $1+1}'`
                        ymax=`awk 'max=="" || $'${column_y}' > max {max=$'${column_y}'} END {printf "%.6f",max}' histo.dat`
                        boundaries=-R${xmin}/${xmax}/0/${ymax}
                        
                        if [[ ${noise_x_ticks} == "" ]] ; then
                                inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.4f",($2-$1)/4}'`
                                inc_x2=`echo "${inc_x}" | awk '{printf "%.4f",($1)/2}'`
                                noise_x_ticks=a${inc_x}f${inc_x2}
                        fi
                        if [[ ${noise_y_ticks} == "" ]] ; then
                                inc_y=`echo "${xmin}" "${xmax}" | awk '{printf "%.4f",($2-$1)/4}'`
                                inc_y2=`echo "${inc_y}" | awk '{printf "%.4f",($1)/2}'`
                                noise_y_ticks=a${inc_y}f${inc_y2}
                        fi
                        
                        bin_width=`awk '{if (NR<=2) print $'${column_x}'}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
                        
                        if [[ ${noise_type} == "3" ]] ; then
                                awk '{print $'${column_x}',$'${column_y}'}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                                                           -B${noise_x_ticks}:"@~l@~":/${noise_y_ticks}:"p(@~l@~|d)": \
                                                                                                           -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                        elif [[ ${noise_type} == "0" ]] ; then
                                awk '{print $'${column_x}',$'${column_y}'}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                                                           -B${noise_x_ticks}:"@~s@~":/${noise_y_ticks}:"p(@~s@~|d)": \
                                                                                                           -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                        fi
                        # Convert ps to png and move file
                        ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
                        
                        echo "...done!"
                        
                fi
        done
        
        if [[ ${noise_type} == "4" ]] ; then
                
                echo "-------------------------------------------------"
                echo "Plotting posterior on noise..."
                
                histo_image=`echo "${postnoise_file}" | awk 'BEGIN {FS="."} ; {print $1}'`
                xmin=`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min}' histo.dat`
                xmax=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max}' histo.dat`
                ymax=`awk 'max=="" || $2 > max {max=$2} END {printf "%.6f",max}' histo.dat`
                boundaries=-R${xmin}/${xmax}/0/${ymax}
                
                if [[ ${noise_x_ticks} == "" ]] ; then
                        inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.4f",($2-$1)/4}'`
                        inc_x2=`echo "${inc_x}" | awk '{printf "%.4f",($1)/2}'`
                        noise_x_ticks=a${inc_x}f${inc_x2}
                fi
                if [[ ${noise_y_ticks} == "" ]] ; then
                        inc_y=`echo "${ymax}" | awk '{printf "%.4f",($1)/4}'`
                        inc_y2=`echo "${inc_y}" | awk '{printf "%.4f",($1)/2}'`
                        noise_y_ticks=a${inc_y}f${inc_y2}
                fi
                
                bin_width=`awk '{if (NR<=2) print $1}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
                awk '{print $1,$2}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                   -B${noise_x_ticks}:"@~s@~":/${noise_y_ticks}:"p(@~s@~|d)": \
                                                                   -Sb${bin_width}u -Ggrey40 -Wthinnest,grey40 -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                # Convert ps to png and move file
                ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
                
                echo "...done!"
                
        fi
fi

# Plot histogram of resistivity posterior
if [[ ${plot_post_value} == "yes" ]] ; then
        
        echo "-------------------------------------------------"
        echo "Plotting posterior on resistivities..."
        
        cp ${file_directory}/${postvalue_file} posterior.dat
        awk '{if (NR<=12) print $0}' ${gmtplot_directory}/${get_posterior}.in | awk -v np=`echo ${number_of_points}` '
                { if(NR==4) print "posterior.dat                     :: Input file";
                else if(NR==8) print "post                           :: Output file (gets renamed to e.g. post_point_1.out)";
                else if(NR==12) print np "                              :: Number of points";
                else if(NR>12) print ;
                else print $0 }
                ' > ${gmtplot_directory}/${get_posterior}1.in
        mv ${gmtplot_directory}/${get_posterior}1.in ${gmtplot_directory}/${get_posterior}.in
        
        for point_number in $(seq 1 1 ${number_of_points}) ; do
                pt_y="p${point_number}_y"
                pt_x="p${point_number}_x"
                echo ${!pt_x} 0 ${!pt_y} >> ${gmtplot_directory}/${get_posterior}.in
        done
        ./${get_posterior}
        
        xmin=${value_min}
        xmax=${value_max}
        
        for point_number in $(seq 1 1 ${number_of_points}) ; do
                echo "   - point ${point_number}"
                histo_image=`echo "${postvalue_file}" | awk 'BEGIN {FS="."} ; {print $1"_'${point_number}'"}'`
                pt_y="p${point_number}_y"
                pt_x="p${point_number}_x"
                cp post_point_${point_number}.out histo.dat
                bin_width=`awk '{if (NR<=2) print $1}' histo.dat | awk 'NR>1{print $1-p} {p=$1}'`
                
                ave=`awk '{print $1*$2}' histo.dat | awk 'BEGIN {FS=" "} ; {sum+=$1} END {print sum*'${bin_width}'}'`
                std=`awk '{print ($1-'${ave}')^2*$2}' histo.dat | awk 'BEGIN {FS=" "} ; {sum+=$1} END {print sqrt(sum*'${bin_width}')}'`
                palette_for_points=stdev_palette.cpt
                if [[ -s ${palette_for_points} ]] ; then
                        color=`awk '{if ($1<='${std}' && $3>'${std}') print $2}' ${palette_for_points}`
                else
                        color=""
                fi
                if [[ ${color} == "" ]] ; then
                        echo "     standard deviation color not found"
                        color=grey40
                fi
                
                ymax=`awk 'max=="" || $2 > max {max=$2} END {print max}' histo.dat`
                boundaries=-R${xmin}/${xmax}/0/${ymax}
                value_y_ticks="p${point_number}_y_ticks"
                value_y_ticks=${!value_y_ticks}
                if [[ ${value_y_ticks} == "" ]] ; then
                        inc_y=`echo "${ymax}" | awk '{printf "%.4f",($1)/4}'`
                        inc_y2=`echo "${inc_y}" | awk '{printf "%.4f",($1)/2}'`
                        value_y_ticks=a${inc_y}f${inc_y2}
                fi
                awk '{print $1,$2}' histo.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                   -B${map_ticks}:"${map_label}":/${value_y_ticks}:"p(@~r@~|d)": \
                                                                   -Sb${bin_width}u -G${color} -Wthinnest,${color} -P -K -Xa${x_corner} -Ya${y_corner} > ${psfile}
                echo "${xmax} ${ymax} 6,1,black 0 RT [${!pt_x} ${!pt_y}]" | ${gmt_prefix} pstext ${boundaries} ${projection_xy} \
                                                                                                  -F+f+a+j -Gwhite -Wblack -N -O -P -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                # Convert ps to png and move file
                ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${histo_image}
        done
        
        echo "...done!"
        
fi

# Plot average number of cells and misfit as a function of iterations
if [[ ${plot_chain_evol} == "yes" ]] ; then
        
        echo "-------------------------------------------------"
        echo "Plotting chain evolution..."
        
        for iter_file in ${ncellsi_file} ${misfiti_file} ; do
                
                iter_image=`echo "${iter_file}" | awk 'BEGIN {FS="."} ; {print $1}'`
                
                if [[ ${iter_file} == "${ncellsi_file}" ]] ; then
                        echo "   - number of cells"
                        awk '{print $1/(10^'${chain_ncell_exp}'),$2}' ${file_directory}/${iter_file} > iter.dat
                        ymin=${ncell_min}
                        ymax=${ncell_max}
                        label="number of cells"
                        chain_exp=${chain_ncell_exp}
                        chain_x_ticks=${chain_ncell_x_ticks}
                        chain_y_ticks=${chain_ncell_y_ticks}
                elif [[ ${iter_file} == "${misfiti_file}" ]] ; then
                        echo "   - misfit"
                        awk '{print $1/(10^'${chain_misfits_exp}'),log($2)/log(10)}' ${file_directory}/${iter_file} > iter.dat
                        ymin=`awk '{if (NR>1) print $0}' iter.dat | awk 'min=="" || $2 < min {min=$2} END {printf "%.6f",min}'`
                        ymax=`awk '{if (NR>1) print $0}' iter.dat | awk 'max=="" || $2 > max {max=$2} END {printf "%.6f",max}'`
                        label="log@-10@-(misfit)"
                        chain_exp=${chain_misfits_exp}
                        chain_x_ticks=${chain_misfits_x_ticks}
                        chain_y_ticks=${chain_misfits_y_ticks}
                fi
                
                xmin=`awk 'min=="" || $1 < min {min=$1} END {printf "%.6f",min}' iter.dat`
                xmax=`awk 'max=="" || $1 > max {max=$1} END {printf "%.6f",max}' iter.dat`
                boundaries=-R${xmin}/${xmax}/${ymin}/${ymax}
                
                if [[ ${chain_x_ticks} == "" ]] ; then
                        inc_x=`echo "${xmin}" "${xmax}" | awk '{printf "%.4f",($2-$1)/4}'`
                        inc_x2=`echo "${inc_x}" | awk '{printf "%.2f",($1)/2}'`
                        chain_x_ticks=a${inc_x}f${inc_x2}
                fi
                if [[ ${chain_y_ticks} == "" ]] ; then
                        inc_y=`echo "${ymin}" "${ymax}" | awk '{printf "%.4f",($2-$1)/4}'`
                        inc_y2=`echo "${inc_y}" | awk '{printf "%.2f",($1)/2}'`
                        chain_y_ticks=a${inc_y}f${inc_y2}
                fi
                
                ${gmt_prefix} psbasemap ${boundaries} ${projection_xy} -B${chain_x_ticks}:"sample (x10@+${chain_exp}@+)":/${chain_y_ticks}:"${label}": \
                                        -K -P -Xa${x_corner} -Ya${y_corner} > ${psfile}
                
                for chain_number in $(seq 1 1 ${number_of_chains}) ; do
                        
                        cd ${inversion_directory}
                        
                        if [[ ${iter_file} == "${ncellsi_file}" ]] ; then
                                
                                ./chain2txt ${ncell_file}${chain_number}
                                
                                awk '{print $1/(10^'${chain_ncell_exp}'),$2}' ${ncell_file}${chain_number}.txt > ${gmtplot_directory}/chain.dat
                                
                                rm ${ncell_file}${chain_number}.txt
                                
                        elif [[ ${iter_file} == "${misfiti_file}" ]] ; then
                                
                                ./chain2txt ${misfits_file}${chain_number}
                                
                                awk '{print $1/(10^'${chain_misfits_exp}'),log($2)/log(10)}' ${misfits_file}${chain_number}.txt > ${gmtplot_directory}/chain.dat
                                
                                rm ${misfits_file}${chain_number}.txt
                                
                        fi
                        
                        cd ${gmtplot_directory}
                        awk '{print $1,$2}' chain.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                           -Wthinnest,blue -P -O -K -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                
                done
                awk '{print $1,$2}' iter.dat | ${gmt_prefix} psxy ${boundaries} ${projection_xy} \
                                                                  -Wthin,red -P -O -Xa${x_corner} -Ya${y_corner} >> ${psfile}
                # Convert ps to png and move file
                ${gmt_prefix} psconvert ${psfile} -A ${plot_format} -F${plot_directory}/${iter_image}
                
        done
        
        echo "...done!"
        
fi

echo "---------------------- DONE ---------------------"

# Output date and time
echo
echo "*************************************************"
echo "Plotting finished:"
date "+%d/%m/%Y %H:%M:%S"
echo "*************************************************"
echo

exit
