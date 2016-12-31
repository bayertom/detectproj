# detectproj

<h2><strong><strong><span style="font-family: Arial;">Software </span></strong><span style="color: #ff0000; font-family: Arial;">detectproj</span><span style="font-family: Arial;">, version 1.1,<br />GNU/GPL projection analysis software for Windows ® 7/8/8.1/10, GNU/Linux and MacOS.</span></strong></h2>
<p>Automated estimation of the map projection and its parameters based on the non-linear optimization...

<ul>
<li>100 map projections are supported,</li>
<li>detection of the projection name and family,</li>
<li>estimation of the normal/transverse/oblique aspect of the projection,</li>
<li>detection of true parallels lat_1, lat_2,</li>
<li>detection of the central parallel shift lon_0,</li>
<li>estimation of the map scale, map rotation (optional),</li>
<li>2 detection methods,</li>
<li>3 optimization techniques,</li>
<li>list of candidate projections sorted by the residuals,</li>
<li>visualization of the detected parameters and residuals,</li>
<li>export reconstructed graticules in DXF,</li>
</ul>


<h3><strong>Running the script</strong></h3>

 1] Basic: set method, test and reference files, detection/optimization methods 

    detectprojv2 +met=nlsm7 test_points.txt reference.txt

    available methods: nlsm7, nlsm8, nmm7, nmm8, dem7, dem8

 		nlsm7		Non-linear least squares optimization, 7 determined parameters (M7)
 		nlsm8		Non-linear least squares optimization, 8 determined parameters (M8), rotation involved
 		nmm7		Nelder-Mead optimization, 7 determined parameters (M7)
 		nmm8		Nelder-Mead optimization, 8 determined parameters (M8), rotation involved
 		dem7		Differential evolution optimization, 7 determined parameters (M7)
 		dem8		Differential evolution optimization, 8 determined parameters (M8), rotation involved

 2] Set method, test and reference files, detection/optimization methods, meridian/parallel increments
 
    Important for DXF with the generated meridians/parallels; setting dlat/dlon increments of meridians/parallels

    detectprojv2 +met=nlsm7 +dlat=10 +dlon=10 test_points.txt reference.txt

 		dlat		Latitude step between two parallels (dlat >=1)
 		dlon		Longitude step between two meridians (dlon>=1)

 3] Set method, test and reference files, detection/optimization methods, meridian/parallel increments, amount of exported graticules to DXF

    detectprojv2 +met=nlsm7 +dlat=10 +dlon=10 +gr=20 test_points.txt reference.txt

 		gr		Amount of best-fit graticules exported to the DXF file (gr<=90)


 Example:

		detectprojv2.exe +met=nlsm7 +dlat=10 +dlon=10 +gr=30 e:\maps\WorldMaps\Seutter\test.txt e:\maps\WorldMaps\Seutter\reference.txt 


<h3>&nbsp;Sample 1: Map of Europe, normal aspect of the projection</h3>
<table style="width: 655px; height: 204px; margin-right: auto; margin-left: auto;">
<tbody>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Title:</span></td>
<td><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;L'Europe sous l'Empire de Charlemagne ou tableau historique de cette partie du monde</span></span></span><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;"> </span> </span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Author: </span></td>
<td><span style="font-size: 12pt; font-family: 'Times New Roman',serif;">&nbsp;Adrien Huber</span></td>
</tr>
<tr>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;Date:</span></td>
<td><span style="font-family: times new roman,times;">&nbsp; <span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">1828,</span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Publisher:&nbsp; </span></td>
<td><span style="font-family: times new roman,times;">&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">A. Brue,</span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Location: </span></td>
<td><span style="font-size: 12pt; font-family: 'Times New Roman',serif;">&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Paris, </span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Type:</span></td>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Atlas Map, </span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Height: </span></td>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;37 cm, </span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Width: </span></td>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">51 cm, </span></td>
</tr>
<tr>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;Scale: </span></td>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">1 : 12,500,000.</span></td>
</tr>
</tbody>
</table>
<p><span style="font-size: 12pt;"><br /><img src="https:web.natur.cuni.cz/~bayertom/images/europe_1_bone_grat_2.jpg" alt="europe 1 bone grat 2" width="885" height="636" style="display: block; margin-left: auto; margin-right: auto;" /></span></p>
<p><span style="font-size: 12pt;">Determined parameters of the projection:<br /></span><span style="font-size: 12pt;"><span style="font-size: 12pt;">Projection: Bonne.<br /></span>Transformed pole position: ?<sub>k</sub> = 90.0<sup>0</sup>, ?<sub>k</sub> = 0.0<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Standard parallels: ?<sub>1</sub> =53. 5<sup>0</sup>, ?2 =53. 5<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Longitude of the central meridian: ?<sub>0</sub> = 14.7<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Arbitrary constant parameter: k = 1.0. </span><br /><span style="font-size: 12pt;"> Auxiliary sphere radius: R'= 0.3m. </span><br /><span style="font-size: 12pt;"> Map scale: S = 20 910 248.</span><br /><span style="font-size: 12pt;"> Angle of rotation: ? = 0.00.</span></p>
<p>&nbsp;</p>
<h3><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Sample 2: World map in hemisphere, transverse aspect of the projection<br /></span></h3>
<table style="width: 655px; height: 204px; margin-right: auto; margin-left: auto;">
<tbody>
<tr>
<td>&nbsp; <span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Title:</span></td>
<td><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Novus Orbis Sive America Meridionalis Et Septentrionalis : divisa per sua regna, provincias et insul, cura et opera</span></span></span></span><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;"> </span> </span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Author: </span></td>
<td><span style="font-size: 12pt; font-family: 'Times New Roman',serif;">&nbsp; <span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Seutter, Matthäus, </span><br /></span></td>
</tr>
<tr>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;Date:</span></td>
<td><span style="font-family: times new roman,times;">&nbsp; <span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">1744</span></span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Publisher:&nbsp; </span></td>
<td><span style="font-family: times new roman,times;">&nbsp; <span style="font-size: 12pt; font-family: 'Times New Roman',serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Seutter, Matthäus, </span></span><br /></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Location: </span></td>
<td><span style="font-size: 12pt; font-family: 'Times New Roman',serif;">&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Augsburg</span>, </span></span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Type:</span></td>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Atlas Map, </span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Height: </span></td>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;58 cm, </span></td>
</tr>
<tr>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">Width: </span></td>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;"><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;50 cm,</span></span></td>
</tr>
<tr>
<td><span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">&nbsp;Scale: </span></td>
<td>&nbsp;<span style="font-size: 12pt; line-height: 107%; font-family: 'Calibri',sans-serif;">1 : 19,000,000</span></td>
</tr>
</tbody>
</table>

<p>&nbsp;<span style="font-family: 'Times New Roman',serif;"> </span><em></em></p>
<p><img src="https:web.natur.cuni.cz/~bayertom/images/america_2_stereo_grat2.jpg" alt="america 2 stereo grat2" width="910" height="785" style="display: block; margin-left: auto; margin-right: auto;" /></p>
<p><span style="font-size: 12pt;"><span style="font-size: 12pt;">Determined parameters of the projection:<br /><span style="font-size: 12pt;">Projection:stereographic.</span></span><br />Transformed pole position: ?<sub>k</sub> = 0.2<sup>0</sup>, ?<sub>k</sub> = -80.3<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Standard parallels: ?<sub>1</sub> =0. 0<sup>0</sup>, ?<sub>2</sub> =0.0<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Longitude of the central meridian: ?<sub>0</sub> = 0.0<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Arbitrary constant parameter: k = 1.0. </span><br /><span style="font-size: 12pt;"> Auxiliary sphere radius: R'= 2.2 m. </span><br /><span style="font-size: 12pt;"> Map scale: S = 39 380 986.</span><br /><span style="font-size: 12pt;"> Angle of rotation: ? = 0.00.</span></p>





=======
"# detectprojv2j" 
>>>>>>> 7fbc3b2cdda953b4c44b50f3080c0e02994e92b8
=======
"# detectprojv2j" 
>>>>>>> 7fbc3b2cdda953b4c44b50f3080c0e02994e92b8
"# detectprojv2j" 
