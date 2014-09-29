
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

<h2>semtools: Methods for Testing Measurement Invariance in Structural Equation Models</h2>

<h3>Testing non-nested structural equation models</h3>

<ul>
  <li><a href="http://arxiv.org/abs/1402.6720">Working paper</a></li>
  <li>Replication materials:<ul>
    <li><a href="sim-vuong.R">Simulation functions</a> for data generation, power evaluation, and power summaries.</li>
    <li><a href="burnout.rda">Burnout data</a> for application.</li>
    <li><a href="replication-vuong.R">Replication script</a> containing code to run and summarize the application and simulations, using the above files.</li>
    <li>Relies on R packages <em>nonnest2</em> to carry out the tests and <em>lavaan</em> for model estimation.
  </ul></li>
</ul>  

<h3>Score-based tests of measurement invariance: Use in practice</h3>

<ul>
  <li>Manuscript published in <a href="http://journal.frontiersin.org/Journal/10.3389/fpsyg.2014.00438/abstract">Frontiers in Psychology 5(438), 1-11</a>.</li>
  <li>Replication materials:<ul>
    <li><a href="mz-frontiers.R">Model estimation functions</a> for simulations.</li>
    <li><a href="sim-frontiers.R">Simulation functions</a> for data generation, power evaluation, and power summaries.</li>
    <li><a href="replication-frontiers.R">Replication script</a> containing code to run and summarize the tutorial and simulations, using the above files. (Note: strucchange 1.5-0 and lavaan 0.5-14 contain code necessary to carry out the tests for general SEMs.)</li>
  </ul></li>
</ul>  

<h3>Testing for measurement invariance with respect to an ordinal variable</h3>

<ul>
  <li>Manuscript accepted for publication in Psychometrika.</li>
  <li><a href="http://econpapers.repec.org/RePEc:inn:wpaper:2012-24">Working paper</a></li>
  <li>Replication materials:<ul>
    <li><a href="estfun-lavaan.R">lavaan extensions</a> containing <tt>estfun()</tt> method for lavaan objects (note: this code has been incorporated into lavaan and is no longer necessary).</li>
    <li><a href="efpFunctional-cat.R">strucchange extensions</a> containing <tt>efpFunctional</tt>s for ordinal measurement invariance tests.</li>
    <li><a href="mz-ordinal.R">Artificial example functions</a> for lavaan model estimation, along with score extraction.</li>
    <li><a href="sim-ordinal.R">Simulation functions</a> for data generation, power evaluation, and power summaries.</li>
    <li><a href="replication-ordinal.R">Replication script</a> containing code to run and summarize the examples and simulations, using the above files.</li>
  </ul></li>
</ul>


<h3>Tests of measurement invariance without subgroups: A generalization of classical methods</h3>

<ul>
  <li>Manuscript published in <a href="http://dx.doi.org/10.1007/S11336-012-9302-4">Psychometrika 78(1), 58-82</a>.</em></li>
  <li><a href="http://www.psychoco.org/2012/slides/Merkle.pdf">Psychoco 2012 presentation</a></li>
  <li>Replication materials:<ul>
    <li><a href="estfun-lavaan.R">lavaan extensions</a> containing <tt>estfun()</tt> method for lavaan objects.</li>
    <li><a href="mz.R">Artificial example functions</a> for lavaan or OpenMx model estimation, along with score extraction.</li>
    <li><a href="sim.R">Simulation functions</a> for data generation, power evaluation, and power summaries.</li>
    <li><a href="replication.R">Replication script</a> containing code to run and summarize the examples and simulations, using the above files.</li>
  </ul></li>
</ul>



<h3>Acknowledgments</h3>

<p> This material is based upon work supported by the U.S. National Science Foundation under Grant No. SES-1061334.  Any opinions, findings and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation (NSF). </p>

<p> The <strong>project summary page</strong> resides <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
