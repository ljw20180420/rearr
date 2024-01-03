

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html><head>


  
  <meta http-equiv="Content-Type" content="text/html; charset=gb2312"><title>kpLogo: k-mer probability logo</title>
  

  
  
  <style>
input[type="submit"]{
/* change these properties to whatever you want */
background-color: #555;
color: #fff;
border-radius: 10px;
font-size:1em;
height:30px;
}
  </style>
  
  <style>
input[type="file"]{
/* change these properties to whatever you want */
background-color: #555;
color: #fff;
border-radius: 10px;
font-size:1em;
}
  </style>
  
  <style>
input[type="text"]{
/* change these properties to whatever you want */
background-color: #aaa;
color: #fff;
}
  </style></head><body style="margin-left: 101px; width: 970px;">

<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" id="title"> <br>
<h1><span class="title"> <i>k</i>pLogo:  <i>k</i>-mer probability logo</span><br> </h1>
</div>

<div style=" border-top: 8px solid DimGray; width: 970px"></div>
<div style=" border-top: 8px solid White; width: 970px" ></div>


</body></html>


<?php // get values of all variables


$foregroundpaste=$_POST['foregroundpaste'];	// foreground sequences
$foregroundfile=$_POST['foregroundfile'];	// forgroundfile
$backgroundpaste=$_POST['backgroundpaste'];	// foreground sequences
$backgroundfile=$_POST['backgroundfile'];	// forgroundfile


$p= scrub_alphanum($_POST['p']);

$small_sample=scrub_alphanum($_POST['small_sample']);

$email=scrub_alphanum($_POST['email']);
$jobname=scrub_alphanum($_POST['jobname']);


$inputtype=scrub_alphanum($_POST['inputtype']); // (none), -ranked, -weighted

$col_seq=scrub_alphanum($_POST['col_seq']);
$col_weight=scrub_alphanum($_POST['col_weight']);

$select_pkmers=scrub_alphanum($_POST['select_pkmers']);
if($select_pkmers != "")
{
	$select_pkmers = " -select ". $select_pkmers;
}
$remove_pkmers=scrub_alphanum($_POST['remove_pkmers']);
if($remove_pkmers != "")
{
    $remove_pkmers = " -remove ". $remove_pkmers;
}


$alphabet=scrub_alphanum($_POST['alphabet']);		// alphabet: DNA, protein, other
$other_alphabet=scrub_alphanum($_POST['other_alphabet']);		// alphabet: DNA, protein, other
if($alphabet == "other")
{
	$alphabet = $other_alphabet;
}

$region_first=scrub_alphanum($_POST['region_first']);
$region_last=scrub_alphanum($_POST['region_last']);


$plottype=scrub_alphanum($_POST['plottype']);

$last_residue=scrub_alphanum($_POST['last_residue']);

$mincount = scrub_alphanum($_POST['mincount']);
if($mincount < 0)
{
	$mincount = 0;
}

$kmer_length = scrub_alphanum($_POST['kmer_length']);

$shift = scrub_alphanum($_POST['shift']);

$degenerate = scrub_alphanum($_POST['degenerate']);
$degenerate_alphabet = scrub_alphanum($_POST['degenerate_alphabet']);

if($degenerate == "other")
{
	if($degenerate_alphabet == "Other: enter all allowed symbols here")
	{
		echo "<font color='red'> Please enter degenerate symbols !</font>";
		exit();	
	}
	$degenerate = " -degenerate $degenerate_alphabet ";
}

if($alphabet != "dna")
{
	$degenerate = "";
}

$background = scrub_alphanum($_POST['background']); // 
$markov_foreground_order = scrub_alphanum($_POST['markov_foreground_order']);
$markov_background_order = scrub_alphanum($_POST['markov_background_order']);
$shuffle_n = scrub_alphanum($_POST['shuffle_n']);
$shuffle_m = scrub_alphanum($_POST['shuffle_m']);
$markov_string = scrub_alphanum($_POST['markov_string']);

if($background == "markov_foreground"){
	$background = " -markov $markov_foreground_order ";
} elseif ($background == "shuffle") {
	$background = " -shuffle $shuffle_n,$shuffle_m ";
} elseif ($background == "bgfile") {
	$background = " -bgfile kpLogo.background.txt ";
} elseif ($background == "markov_background") {
	$background = " -markov $markov_background_order -bgfile kpLogo.background.txt ";
} elseif ($background == "markov_string") 
{
	$background = " -markov $markov_string";
}

$pseudo = scrub_alphanum($_POST['pseudo']);

$maxFrac = scrub_alphanum($_POST['maxFrac']);

	
$startPos = scrub_alphanum($_POST['startPos']);
$colorblind = scrub_alphanum($_POST['colorblind']);

$stack_order = scrub_alphanum($_POST['stack_order']);


// clean up folders older than 3 days
$clean = "find ./files_to_be_removed_10_days_after_creation -mtime +10 -exec rm -rf {} \;";
exec('nohup '. $clean . ' > /dev/null 2>&1 &');


// make a folder with randome name for each job
$jobID = substr(str_shuffle(md5(time())),0,20);
 
$tmpfolder = "./files_to_be_removed_10_days_after_creation/".$jobID;

$oldmask = umask(0);
if (!mkdir($tmpfolder, 0777, true)) {
    die('Failed to create folders...');
}

// enter this folder
chdir($tmpfolder);

// input
if (strlen($foregroundpaste) > 0) // sequence pasted, save to file
{
	$h = fopen("kpLogo.input.txt", 'w');
	fwrite($h, $foregroundpaste);
	fclose($h); 
}
else
{
	$file1 = $_FILES['foregroundfile'];
	if(move_uploaded_file($file1['tmp_name'], "kpLogo.input.txt")) 
	{
		//echo "The file ".basename( $_FILES['uploadedfile']['name'])." has been uploaded";
	}
	else
	{
		echo "<font color='red'>Please paste or upload your input!</font>";
		exit();
	}
}

// background
if ($background == " -bgfile kpLogo.background.txt " || $background == " -markov $markov_background_order -bgfile kpLogo.background.txt ") {
	if (strlen($backgroundpaste) > 0) // sequence pasted, save to file
	{
		$hb = fopen("kpLogo.background.txt", 'w');
		fwrite($hb, $backgroundpaste);
		fclose($hb); 
	}
	else 
	{
		$file2 = $_FILES['backgroundfile'];
		if(move_uploaded_file($file2['tmp_name'], "kpLogo.background.txt"))
		{
		}
		else 
		{
			echo "<font color='red'>Fail to upload background file!</font>";
			exit();
		}
	}
}


//

$command = "kpLogo kpLogo.input.txt -o kpLogo.output $inputtype -seq $col_seq -weight $col_weight -alphabet $alphabet  $kmer_length $shift $background -startPos $startPos $degenerate $colorblind -minCount $mincount -pseudo $pseudo -region $region_first,$region_last $select_pkmers $remove_pkmers $plottype $small_sample -pc $p $last_residue $stack_order -fix $maxFrac";

file_put_contents("../../num.job.each.day", date("Y-m-d")."\n", FILE_APPEND | LOCK_EX);
$totalLines=intval(exec('wc -l ../../num.job.each.day'));
$jobID="kpLogo-".$totalLines;


//email
//$subject = "kpLogo results available: $jobID ($jobname)";
//$url = str_replace("submit.php","$tmpfolder/result.php","http://$_SERVER[HTTP_HOST]$_SERVER[REQUEST_URI]");
//$content = "Your kpLogo job $jobID ($jobname) is finished and results are available here for *** 10 days ***: \r\n\r\n $url";

//$result = exec('nohup ../../'. $command . ' -email '. $email . ' -subject "'. $subject. '" -content "'. $content. '" >> log 2>&1 &');
$result = exec('nohup ../../'. $command .' >> log 2>&1 &');

umask($oldmask);

exec("cp ../../result.php ./");


// save relevent information
file_put_contents("job.info", $folder."\n", FILE_APPEND | LOCK_EX);
file_put_contents("job.info", $jobname."\n", FILE_APPEND | LOCK_EX);
file_put_contents("job.info", $command."\n", FILE_APPEND | LOCK_EX);
file_put_contents("job.info", $jobID."\n", FILE_APPEND | LOCK_EX);

file_put_contents("submission_notification_not_sent", $email."\n", FILE_APPEND | LOCK_EX);
file_put_contents("submission_notification_not_sent", "This file will be deleted once the job is finished and the email is sent\n", FILE_APPEND | LOCK_EX);




header("Location: $tmpfolder/result.php");


function scrub_alphanum($str){
  $str = preg_replace("/\.+/",".",$str); # All repeating dots to single dots
  $str = preg_replace("/[^@0-9a-zA-Z-_:,. ]/","",$str); # Remove all non alphanum, comma, period, at, hyphen
  return $str;
}


?>



