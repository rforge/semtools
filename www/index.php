
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://R-Forge.R-project.org/"><img src="<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<h2>semtools: Methods for Testing Measurement Invariance in Structural Equation Models</h1>

<h3>Tests of measurement invariance without subgroups: A generalization of classical methods</h2>

<ul>
  <li><a href="http://econpapers.repec.org/RePEc:inn:wpaper:2011-09">Working paper</a>, revised version accepted for publication in <em>Psychometrika</em></li>
  <li><a href="http://www.psychoco.org/2012/slides/Merkle.pdf">Psychoco 2012 presentation</a></li>
  <li>Replication materials:<ul>
    <li> <a href="mz.R">Artificial example functions</a> for lavaan or OpenMx model estimation, along with score extraction.</li>
    <li> <a href="sim.R">Simulation functions</a> for data generation, power evaluation, and power summaries.</li>
    <li> <a href="replication.R">Replication script</a> containing code to run and summarize the examples and simulations, using the above files.</li>
  </ul></li>
</ul>

<h3>Testing for measurement invariance with respect to an ordinal variable</h2>

<ul>
  <li><a href="http://econpapers.repec.org/RePEc:inn:wpaper:2012-24">Working paper</a>, forthcoming</li>
</ul>


<h2>Acknowledgments</h2>

<p> This material is based upon work supported by the U.S. National Science Foundation under Grant No. SES-1061334.  Any opinions, findings and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation (NSF). </p>

<p> The <strong>project summary page</strong> resides <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
