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
<h2>&nbsp;</h2>
<h2>Command-line version of detectproj in C++</h2>
<p>The generic C++ version&nbsp;(C++11 support) of the software detectprojv2j is available from git repository. <br />The basic parameters can be set using the command-line.<br />All methods and optimization techniques are involved.<br />Supported compilers: g++, msvc2015.</p>
<h3>Basic usage</h3>
<p>Set method, test and reference files, detection/optimization methods. Available methods: nlsm7, nlsm8, nmm7, nmm8, dem7, dem8:</p>
<table>
<tbody>
<tr>
<td><strong>Method</strong></td>
<td><strong>Amount of parameters</strong></td>
<td><strong>Optimization</strong></td>
<td><strong>Description</strong></td>
</tr>
<tr>
<td><em>nlsm7</em></td>
<td style="text-align: center;">7</td>
<td style="text-align: center;">Local</td>
<td>Non-linear least squares optimization based on hybrid BFGS.</td>
</tr>
<tr>
<td><em>nlsm8</em></td>
<td style="text-align: center;">8</td>
<td style="text-align: center;">Local</td>
<td>Non-linear least squares optimization based on hybrid BFGS, map rotation involved.</td>
</tr>
<tr>
<td><em>nmm7</em></td>
<td style="text-align: center;">7</td>
<td style="text-align: center;">Global</td>
<td>Nelder-Mead optimization (Simplex method).</td>
</tr>
<tr>
<td><em>nmm8</em></td>
<td style="text-align: center;">&nbsp;8</td>
<td style="text-align: center;">Global</td>
<td>Nelder-Mead optimization (Simplex method), map rotation involved.</td>
</tr>
<tr>
<td><em>dem7</em></td>
<td style="text-align: center;">7</td>
<td style="text-align: center;">Global</td>
<td>Differential evolution optimization.</td>
</tr>
<tr>
<td><em>dem8</em></td>
<td style="text-align: center;">8</td>
<td style="text-align: center;">Global</td>
<td>Differential evolution optimization, map rotation involved.</td>
</tr>
</tbody>
</table>
<p>Example:</p>
<p><span style="font-family: courier new,courier;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; detectprojv2 +met=nlsm7 e:\maps\Seutter\test.txt e:\maps\Seutter\reference.txt<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; detectprojv2 +met=nmm7 e:\maps\Seutter\test.txt e:\maps\Seutter\reference.txt<br /></span></p>
<h3>Set meridian/parallel increment</h3>
<p>Set method, test and reference files, detection/optimization methods, meridian/parallel increments.&nbsp; This option is important for the DXF files with the generated graticule for setting dlat/dlon increments of meridians/parallels.</p>
<table>
<tbody>
<tr>
<td><strong>Parameter</strong></td>
<td><strong>Default</strong></td>
<td><strong>Description</strong></td>
</tr>
<tr>
<td><em>dlat</em></td>
<td style="text-align: center;">10</td>
<td>Latitude step between two adjacent parallels</td>
</tr>
<tr>
<td><em>dlon</em></td>
<td style="text-align: center;">10</td>
<td>Longitude step between two adjacent meridians</td>
</tr>
</tbody>
</table>
<p>Example:</p>
<p><span style="font-family: courier new,courier;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; detectprojv2 +met=nlsm7 +dlat=15 +dlon=30 e:\maps\Seutter\test.txt e:\maps\Seutter\reference.txt<br /></span></p>
<h3>Set amount of exported graticules</h3>
<p>Set method, test and reference files, detection/optimization methods, meridian/parallel increments, amount of generated graticules exported to the DXF format.</p>
<table>
<tbody>
<tr>
<td><strong>Parameter</strong></td>
<td><strong>Default</strong></td>
<td><strong>Desciption </strong></td>
</tr>
<tr>
<td><em>gr</em></td>
<td style="text-align: center;">20</td>
<td>Amount of best-fit graticules exported to the DXF file; gr&lt;=90.</td>
</tr>
</tbody>
</table>
<p>&nbsp;Example:</p>
<p><span style="font-family: courier new,courier;">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; detectprojv2 +met=nlsm7 +dlat=15 +dlon=30&nbsp;+gr=30 e:\maps\Seutter\test.txt e:\maps\Seutter\reference.txt<br /></span></p>




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
<p><span style="font-size: 12pt;"><br /><img src="https://web.natur.cuni.cz/~bayertom/images/europe_1_bone_grat_2.jpg" alt="europe 1 bone grat 2" width="885" height="636" style="display: block; margin-left: auto; margin-right: auto;" /></span></p>
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
<p><img src="https://web.natur.cuni.cz/~bayertom/images/america_2_stereo_grat2.jpg" alt="america 2 stereo grat2" width="910" height="785" style="display: block; margin-left: auto; margin-right: auto;" /></p>
<p><span style="font-size: 12pt;"><span style="font-size: 12pt;">Determined parameters of the projection:<br /><span style="font-size: 12pt;">Projection:stereographic.</span></span><br />Transformed pole position: ?<sub>k</sub> = 0.2<sup>0</sup>, ?<sub>k</sub> = -80.3<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Standard parallels: ?<sub>1</sub> =0. 0<sup>0</sup>, ?<sub>2</sub> =0.0<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Longitude of the central meridian: ?<sub>0</sub> = 0.0<sup>0</sup>.</span><br /><span style="font-size: 12pt;"> Arbitrary constant parameter: k = 1.0. </span><br /><span style="font-size: 12pt;"> Auxiliary sphere radius: R'= 2.2 m. </span><br /><span style="font-size: 12pt;"> Map scale: S = 39 380 986.</span><br /><span style="font-size: 12pt;"> Angle of rotation: ? = 0.00.</span></p>





=======
"# detectprojv2j" 
>>>>>>> 7fbc3b2cdda953b4c44b50f3080c0e02994e92b8
=======
"# detectprojv2j" 
>>>>>>> 7fbc3b2cdda953b4c44b50f3080c0e02994e92b8
"# detectprojv2j" 
