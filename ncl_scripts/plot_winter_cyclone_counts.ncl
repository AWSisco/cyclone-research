begin

;The data is structured with 25 rows and 144 columns
nrows = 25
ncols= 144
data = asciiread("../subset_cyclones/winter_cyclone_counts.txt", (/nrows, ncols/), "integer")

lat = fspan(20,80,nrows)
lon = fspan(0,357.5,ncols)

data!0 = "lat"
data!1 = "lon"
data&lat = lat
data&lon = lon
data&lat@units = "degrees_north"
data&lon@units = "degrees_east"

data = where(data.eq.0,data@_FillValue,data)

pac_data = data({20:60},{130:240})
atl_data = lonFlip(data)
atl_data := atl_data({20:60},{-100:10})

res = True
res@cnFillOn = True
res@cnLinesOn = False
res@cnFillMode = "RasterFill"
res@gsnAddCyclic = False
res@gsnMaximize = True
res@mpDataBaseVersion = "MediumRes"
res@gsnDraw = False
res@gsnFrame = False
res@mpGridAndLimbOn = True
res@mpGridLonSpacingF = 10
res@mpGridLatSpacingF = 5
res@mpGridLineColor = "grey"
res@mpGridLineDashPattern = 5
res@gsnMinorLatSpacing = 5
res@tmYROn = False ;No right tickmarks
res@tmXTOn = False ;No top tickmarks

res_atl = res

atl_colors = read_colormap_file("MPL_Blues") ;Get all colors
atl_inds = span_color_indexes("MPL_Blues", 10) ;Get indices for 10 colors
atl_colors := atl_colors(atl_inds(1:9)-2, :) ;Extract the 9 colors I want to use
transparent = new((/1,4/), "float") ;Create empty array
transparent(0,:) =  (/255, 255, 255, 0/) / 255.0 ;Fill empty array with transparent RGB
res_atl@cnFillColors = array_append_record(transparent, atl_colors, 0) ;Build the fill color array

res_atl@cnLevelSelectionMode = "ExplicitLevels"
res_atl@cnLevels = (/1,4,8,12,16,20,24,28,32/)

res_atl@mpMinLatF = 20
res_atl@mpMaxLatF = 60
res_atl@mpMinLonF = 260
res_atl@mpMaxLonF = 370
res_atl@mpCenterLonF = 315
res_atl@gsnLeftString = "ATL"
res_atl@gsnRightString = "n = 601"

atl_center_ypts_1 = (/40,40,45,45,40/)
atl_center_xpts_1 = (/-75,-15,-15,-75,-75/)

;================================================

res_pac = res

pac_colors = read_colormap_file("MPL_Blues") ;Get all colors
pac_inds = span_color_indexes("MPL_Blues", 10) ;Get indices for 10 colors
pac_colors := pac_colors(pac_inds(1:9)-2, :) ;Extract the 9 colors I want to use
transparent = new((/1,4/), "float") ;Create empty array
transparent(0,:) =  (/255, 255, 255, 0/) / 255.0 ;Fill empty array with transparent RGB
res_pac@cnFillColors = array_append_record(transparent, pac_colors, 0) ;Build the fill color array

res_pac@cnLevelSelectionMode = "ExplicitLevels"
res_pac@cnLevels = (/1,5,10,15,20,25,30,35,40/)

res_pac@mpMinLatF = 20
res_pac@mpMaxLatF = 60
res_pac@mpMinLonF = 130
res_pac@mpMaxLonF = 240
res_pac@mpCenterLonF = 190
res_pac@gsnLeftString = "WNP"
res_pac@gsnRightString = "n = 738"

pac_west_ypts_0 = (/40,40,45,45,40/)
pac_west_xpts_0 = (/135,195,195,135,135/)
redBoxRes = True
redBoxRes@gsLineThicknessF = 6
redBoxRes@gsLineColor = "red"

wks = gsn_open_wks("x11", "cyclone_counts")
plot = new(2,graphic)
plot(0) = gsn_csm_contour_map(wks, pac_data, res_pac)
pac_west_box_0 = gsn_add_polyline(wks, plot(0), pac_west_xpts_0, pac_west_ypts_0, redBoxRes)
plot(1) = gsn_csm_contour_map(wks, atl_data, res_atl)
atl_center_box_1 = gsn_add_polyline(wks, plot(1), atl_center_xpts_1, atl_center_ypts_1, redBoxRes)

gsn_panel(wks,plot,(/2,1/),False)

end
