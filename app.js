var inputSeq = "";
$("#resultsHead").hide();
$(".references").hide();
$("#ecInput").hide();
$("#reload").hide();
$("#subSeq").on("click", function(){
inputSeq = $("textarea:input[name=inputSeq]").val();
inputSeqNoReturn = inputSeq.replace(/(\r\n|\n|\r)/gm,"");
checkInput(inputSeqNoReturn);
});

//========Check if there is an input========
function checkInput(inputSeqNoReturn){
  var seqLength = inputSeqNoReturn.length;
  var sequence = inputSeqNoReturn.toUpperCase();

  if (seqLength === 0){
  $(".introSection").remove();
  $("#errorMessage").append("No sequence entered. <br>Please enter a DNA or Protein sequence.");

  $("#subSeq").on("click", function(){
  $("#errorMessage").remove();
  $(".errMessage").append("<div id= 'errorMessage'></div>");
  $("#errorMessage").append("No sequence entered.<br>Please enter a DNA or Protein sequence.");
  })} else {
      $(".introSection").hide();
      checkInputType(sequence, seqLength);};
      };

//========Check sequence input types========
function checkInputType(sequence, seqLength){
  if ($("#proteinButton").is(":checked")){
    $(".frontPage").toggle();
    aminoComp(sequence, seqLength);
  } else if ($("#dnaButton").is(":checked")){
    $(".frontPage").toggle();
    dnaComp(sequence, seqLength);
    }  else {
        $(".introSection").remove();
        $("#errorMessage").append("<h2>You did not select a sequence type</h2>");
        $("#subSeq").on("click", function(){
        $("#errorMessage").remove();
        $(".errMessage").append("<div id= 'errorMessage'></div>");
        $("#errorMessage").append("You did not select a sequence type");
        });
      }
};

//=======Protein Module for Computing paramenters========
function aminoComp(sequence, seqLength){
  var mw = 18.015;
  var monMW = 18.010565;
  var molPercent = 0;
  var hydrop = 0;
  var atoms = 0
  var nonAAres = [];
  var nonAApos = [];

var aminoacids = {
  A : ["Alanine", 71.0788, 0, 1.800, 71.03711, 13],
  R : ["Arginine", 156.1875, 0, -4.500, 156.10112, 26],
  N : ["Asparagine", 114.1038, 0, -3.500, 114.04293, 17],
  D : ["Aspartic acid", 115.0886, 0, -3.500, 115.02695, 16],
  C : ["Cysteine", 103.1388, 0, 2.500, 103.00919, 14],
  Q : ["Glutamine", 128.1307, 0, -3.500, 128.05858, 20],
  E : ["Glutamic acid", 129.1155, 0, -3.500, 129.04260, 19],
  G : ["Glycine", 57.0519, 0, -0.400, 57.02147, 10],
  H : ["Histidine", 137.1411, 0, -3.200, 137.05891, 20],
  I : ["Isoleucine", 113.1594, 0, 4.500, 113.08407, 22],
  L : ["Leucine", 113.1594, 0, 3.800, 113.08407, 22],
  K : ["Lysine", 128.1741, 0, -3.900, 128.09497, 24],
  M : ["Methionine", 131.1926, 0, 1.900, 131.04049, 20],
  F : ["Phenylalanine", 147.1766, 0, 2.800, 147.06842, 23],
  P : ["Proline", 97.1167, 0, -1.600, 97.05277, 17],
  S : ["Serine", 87.0782, 0, -0.800, 87.03203, 14],
  T : ["Threonine", 101.1051, 0, -0.700, 101.04768, 17],
  W : ["Tryptophan", 186.2132, 0, -0.900, 186.07932, 27],
  Y : ["Tyrosine", 163.1760, 0, -1.300, 163.06333, 24],
  V : ["Valine", 99.1326, 0, 4.200, 99.06842, 19]
};

for (var i = 0; i< seqLength; i++){
  var res = sequence[i];
    if (!aminoacids.hasOwnProperty(res)){
      nonAAres.push(res);
      nonAApos.push(i+1);
    } else {
        for (var aa in aminoacids){
            if (aa === res){
              mw += aminoacids[aa][1];
              aminoacids[aa][2] += 1;
              hydrop += aminoacids[aa][3];
              monMW += aminoacids[aa][4];
              atoms += aminoacids[aa][5];
            };
        };
      };
};

// % aliphatic index calculation
var percentA = 100*(aminoacids.A[2]/seqLength);
var percentV = 100*(aminoacids.V[2]/seqLength);
var percentI = 100*(aminoacids.I[2]/seqLength);
var percentL = 100*(aminoacids.L[2]/seqLength);
var aliphIndex = percentA + (2.9 * percentV) + (3.9 * (percentI + percentL));

// GRAVY calculation
var gravy = hydrop / seqLength;

// Function to deal with non-standard amino acids
if (nonAAres.length > 0){
      nonAA (nonAAres,nonAApos);
};

// Section to break sequence into blocks
if (nonAAres.length === 0){
  if ($("#proteinButton").is(":checked")){
      $("#resultsHead").show().append("<h3>RESULTS OF YOUR SEQUENCE ANALYSIS</h3>");
      $("#proteinSeqSubHead span").append("<br>PROTEIN ANALYSIS")
      $("#proteinSeqSubHead").append("<br><strong>Sequence</strong><br>The submitted protein sequence has " + seqLength + " amino acids:<br>");
  };
  for (i = 1; i <= seqLength; i++){
    if (i %10 > 0){
      $("#sequence").append(sequence[i-1])
    } else
    $("#sequence").append(sequence[i-1] + "     ")
      if (i % 50 === 0) {
        $("#sequence").append("   " + i + "<br>")
      };
        if (i=== seqLength && i % 50 != 0){
        $("#sequence").append("   " + seqLength)
        }
  };

// Output of protein molecular weights
$("#mw").append("<br><strong>Molecular Weight</strong>")
$("#mw").append("<br>The average molecular weight is: " + mw.toFixed(4));
$("#mw").append("<br>The monisotopic molecular weight is: " + monMW.toFixed(4));
$("#mw").append("<br>The total number of atoms is: " + (atoms - ((seqLength-1)*3)));

$("#composition").append("<br><strong>Amino acid composition</strong><br>");
$("#aaTable thead").append("<th>Amino Acid</th><th>Count</th><th>% Total</th>");

// Output of protein composition table
for (i in aminoacids){
$("#aaTable tbody").append("<tr><td style='text-align:right'>"+ aminoacids[i][0] + " (" + i + ")" + "</td><td>"+ aminoacids[i][2] + "</td><td>" +((aminoacids[i][2]/seqLength)*100).toFixed(2) + "</td><</tr>")
};

// Net charge calculation
var negCharge = (aminoacids.D[2] + aminoacids.E[2]);
var posCharge = (aminoacids.R[2] + aminoacids.K[2] + aminoacids.H[2]);
var hydrophobic = (aminoacids.A[2] + aminoacids.G[2] + aminoacids.I[2] + aminoacids.L[2] + aminoacids.M[2] + aminoacids.F[2] + aminoacids.P[2] + aminoacids.W[2] + aminoacids.V[2]);
var polar = (aminoacids.N[2]+aminoacids.C[2]+aminoacids.Q[2]+aminoacids.S[2]+aminoacids.T[2]+aminoacids.Y[2]);

// Residue character composition calculation and output
$("#chargeHydro").append("<strong>Charge and hydrophobicity<br></strong>");
$("#chargeHydro").append("Negatively charged residues (D, E) = " + negCharge +  " (" + ((negCharge/seqLength)*100).toFixed(2) + "%)<br>");
$("#chargeHydro").append("Positively charged residues (K, R, H) = " + posCharge +  " (" + ((posCharge/seqLength)*100).toFixed(2) + "%)<br>");
$("#chargeHydro").append("Polar residues (C, S, Q, N, T, Y) = " + polar + " (" + ((polar/seqLength)*100).toFixed(2) + "%)<br>");
$("#chargeHydro").append("Hydrophobic residues (A, G, I, L, M, P, F, W, V) = " + hydrophobic + " (" + ((hydrophobic/seqLength)*100).toFixed(2) + "%)<br>");

// Hydrophocity parameters calculation
$("#chargeHydro").append("<br>The aliphatic index is: " +  aliphIndex.toFixed(2));
$("#chargeHydro").append("<br>The Grand Average of Hydropathicity (GRAVY) is: " +  gravy.toFixed(3));
exCoef(aminoacids, mw)
pICalc(aminoacids)}
};

// Function to handle non-standard amino acid codes
function nonAA (nonAAres,nonAApos){
  $("#notNormal span").append("<h2>WARNING!</h2>");
  if (nonAAres.length == 1){
  $("#notNormal").append("Your sequence contains the following non-standard amino acid code: " + nonAAres.toString());
  $("#notNormal").append("<br>It is located at position: " + nonAApos.toString());
  } else {
    $("#notNormal").append("Your sequence contains the following non-standard amino acid codes: " + nonAAres.toString());
    $("#notNormal").append("<br>They are located at positions: " + nonAApos.toString());
    };
  $("#notNormal").append('<br>If you want to enter a revised sequence, click "Enter"');
  $("#notNormal").append("<br><button type='button' id = 'reSub' name='reSub' class='btn btn-danger'>RESUBMIT</button>");

  $("#reSub").on("click", function(){
    $("#notNormal").remove();
    $(".notNorm").append("<div id='notNormal'></div>");

    $("#errorMessage").remove();
    $(".errMessage").append("<div id= 'errorMessage'></div>");

    $(".frontPage").toggle()});
    $("#textInput").val("");
};

// pI calculation
// pI using Wikipedia pKa values
function pICalc(aminoacids){
  var pH = 6.5;
  var pHprev = 0.0;
  var pHnext = 14.0;
  var temp = 0.0;
  var bpH = 6.5;
  var bpHprev = 0.0;
  var bpHnext = 14.0;
  var btemp = 0.0;

do {
  var QN1 = -1/(1+Math.pow(10,(3.65-pH)));
  var QN2 = -1*(aminoacids.D[2])/(1+Math.pow(10,(3.9-pH)));
  var QN3 = -1*(aminoacids.E[2])/(1+Math.pow(10,(4.07-pH)));
  var QN4 = -1*(aminoacids.C[2])/(1+Math.pow(10,(8.18-pH)));
  var QN5 = -1*(aminoacids.Y[2])/(1+Math.pow(10,(10.46-pH)));
  var QP1 = (aminoacids.H[2])/(1+Math.pow(10,(pH-6.04)));
  var QP2 = 1/(1+Math.pow(10,(pH-8.2)));
  var QP3 = (aminoacids.K[2])/(1+Math.pow(10,(pH-10.54)));
  var QP4 = (aminoacids.R[2])/(1+Math.pow(10,(pH-12.48)));
  var NQ = QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4;

  if (NQ < 0){
    temp = pH;
    pH = pH - ((pHnext-pH)/2);
    pHnext = temp
  } else{
    temp = pH;
    pH = pH + ((pHnext-pH)/2);
    pHprev = temp;
    };
}

while (pH-pHprev > 0.01 && pHnext-pH > 0.01);

// pI using Blenquist pKa values
do {
  var bQN1 = -1/(1+Math.pow(10,(3.55-pH)));
  var bQN2 = -1*(aminoacids.D[2])/(1+Math.pow(10,(4.05-bpH)));
  var bQN3 = -1*(aminoacids.E[2])/(1+Math.pow(10,(4.45-bpH)));
  var bQN4 = -1*(aminoacids.C[2])/(1+Math.pow(10,(9.0-bpH)));
  var bQN5 = -1*(aminoacids.Y[2])/(1+Math.pow(10,(10.0-bpH)));
  var bQP1 = (aminoacids.H[2])/(1+Math.pow(10,(bpH-5.98)));
  var bQP2 = 1/(1+Math.pow(10,(bpH-7.5)));
  var bQP3 = (aminoacids.K[2])/(1+Math.pow(10,(bpH-10.0)));
  var bQP4 = (aminoacids.R[2])/(1+Math.pow(10,(bpH-12.0)));
  var bNQ = bQN1+bQN2+bQN3+bQN4+bQN5+bQP1+bQP2+bQP3+bQP4;

    if (bNQ < 0){
      btemp = bpH;
      bpH = bpH - ((bpHnext-bpH)/2);
      bpHnext = btemp} else{
      btemp = bpH;
      bpH = bpH + ((bpHnext-bpH)/2);
      bpHprev = btemp;};
}
while ( bpH-bpHprev > 0.01 && bpHnext-bpH > 0.01 );

$("#pI").append("<strong>Isoelectric Point</strong><br>");
$("#pI").append(bpH.toFixed(2) + " using amino acid pKa values from Bjellqvist et al<br>");
$("#pI").append(pH.toFixed(2) + " using amino acid pKa values from Wikipedia");
};

// Extinction coefficient calculation
function exCoef(aminoacids, mw){
  var trp = aminoacids.W[2];
  var tyr = aminoacids.Y[2];
  if (aminoacids.C[2] >= 2){
    var cys = Math.floor(((aminoacids.C[2])/2))}
    else {var cys = 0};

  var excoefNoCys = (trp*5500) + (tyr*1490);
  var excoefCys = (trp*5500) + (tyr*1490) + (cys*125);
  var concMolar = 1/mw
  var mgmlOdNoCys = (excoefNoCys*concMolar);
  var mgmlOdCys = (excoefCys*concMolar);

  $("#ec").append("<strong>Extinction Coefficient</strong><br>");
  $("#ec").append("The extinction coefficient is: " + excoefNoCys);
  $("#ec").append("<br>A 1mg/ml (" + (concMolar/Math.pow(10,-6)).toFixed(3) + " uM) solution of your protein has an A280nm of: " + mgmlOdNoCys.toFixed(2) + "<br>");

  if (aminoacids.C[2] >= 2){
    $("#ec").append("<br>If all cysteines are disulphide bonded, the extinction coefficient is: " + excoefCys);
    $("#ec").append("<br>A 1mg/ml (" + (concMolar/Math.pow(10,-6)).toFixed(3) + " uM) solution of your protein has an A280nm of: " + mgmlOdCys.toFixed(2) + "<br>");
  };
  $("#ecInput").show();
  $("#reload").show();
  $('#A280button').on('click', function(event){
    event.preventDefault();
    var absorb = parseFloat($('#AbsorbInput').val());

    var inputConcMgMlNoCys = (absorb / mgmlOdNoCys).toFixed(4);
    var inputConcMolarNoCys = ((inputConcMgMlNoCys / mw) * Math.pow(10,6)).toFixed(2);

    if (excoefNoCys === 0 || excoefCys === 0){
    $("#odResult").append("It is not possible to determine the concentration of your protein from its A280.<br>It has an extinction coefficient of 0 (i.e. no tryptophans, tyrosines or pairs of cysteines).")}
    else {
    $("#odResult").append("The concentration of your protein with A280 = " + absorb + " is:<br>" + inputConcMgMlNoCys + " mg/ml  (" +inputConcMolarNoCys + " uM)<br>");
    if (aminoacids.C[2] >= 2){
      var inputConcMgMlCys = (absorb / mgmlOdCys).toFixed(4);
      var inputConcMolarCys = ((inputConcMgMlCys / mw) * Math.pow(10,6)).toFixed(2);
      $("#odResult").append("If the cyteines are all disulphide bonded, the concentration is:<br>" + inputConcMgMlCys + " mg/ml  (" +inputConcMolarCys + " uM)<br><br>");
    };
    };
    $("#AbsorbInput").val("");
    });

    $('#reload').on('click', function() {
       location.reload();
    });
    $(".references").show();
};


// =======DNA analysis module=============
function dnaComp(sequence, seqLength){
  var nonDNAbase = [];
  var nonDNApos = [];
  var baseList = ["A", "T", "C", "G"];
  for (var i=0; i < seqLength; i++){
    var base = sequence[i];
    if (!baseList.includes(base)){
      nonDNAbase.push(base);
      nonDNApos.push(i+1);
    };
  };
  if (nonDNAbase.length > 0){
    nonDNA (nonDNAbase,nonDNApos);
  } else {
  dnaAnal(sequence, seqLength);
};
};
// Function to handle non-standard bases
function nonDNA (nonDNAbase,nonDNApos){
  $("#notNormal span").append("<h2>WARNING!</h2>");
  if (nonDNAbase.length == 1){
  $("#notNormal").append("Your sequence contains the following character that is not a base: " + nonDNAbase.toString() + ". ");
  $("#notNormal").append("<br>It is located at position: " + nonDNApos.toString() + ".");
  } else {
  $("#notNormal").append("Your sequence contains the following characters that are not bases: " + nonDNAbase.toString() + ". ");
  $("#notNormal").append("<br>They are located at positions: " + nonDNApos.toString() + ".");
    };

  $("#notNormal").append('<br>If you want to enter a revised sequence, click "RESUBMIT:"');
  $("#notNormal").append("<br><button type='button' id = 'reSub' name='reSub' class='btn btn-danger'>RESUBMIT</button>");

  $("#reSub").on("click", function(){
  $("#notNormal").remove();
  $(".notNorm").append("<div id='notNormal'></div>");

  $("#errorMessage").remove();
  $(".errMessage").append("<div id= 'errorMessage'></div>");

  $(".frontPage").toggle()});
  $("#textInput").val("");
};

// Function to analyse DNA sequence composition
function dnaAnal(sequence, seqLength){
  var stop = [];
  var codons = [];
  var gComp, cComp, aComp, tComp;
  gComp = cComp = aComp = tComp = 0;

  for (var i=0; i < seqLength; i++){
    if (sequence[i] === "G"){gComp++}
    else if (sequence[i] === "C") {cComp++}
    else if (sequence[i] === "A"){aComp++}
    else {tComp++}
  };

  for (var i=0; i < seqLength; i=i+3){
    var codon = sequence.slice(i,i+3);
    codons.push(codon);
    if (codon === "TGA" || codon==="TAG" || codon === "TAA"){
    stop.push(codons.length);
    };
  };

  $("#resultsHead").show().append("<h2>RESULTS OF YOUR SEQUENCE ANALYSIS</h2>");
  $("#dnaSeqSubHead span").append("<br><strong>DNA ANALYSIS</strong><br>");
  $("#dnaSeqSubHead").append("The DNA sequence entered has " + seqLength + " bases: <br>");

  for (var i = 1; i < (codons.length+1); i++){
    $("#dnaSequence").append(codons[i-1] + "     ");
    if (i%20 === 0){
      $("#dnaSequence").append("    " + (i*3) + "<br>");
    };
  };

  $("#dnaSequence").append("    " + (codons.length*3) + "<br>");

  var gcRatio = ((gComp+cComp) / (seqLength))*100;
  $("#dnaSeqComp").append("<br>% G: " + ((gComp/seqLength)*100).toFixed(2) + "<br>");
  $("#dnaSeqComp").append("% C: " + ((cComp/seqLength)*100).toFixed(2) + "<br>");
  $("#dnaSeqComp").append("% A: " + ((aComp/seqLength)*100).toFixed(2) + "<br>");
  $("#dnaSeqComp").append("% T: " + ((tComp/seqLength)*100).toFixed(2) + "<br>");
  $("#dnaSeqComp").append("The %GC:%AT ratio is: " + gcRatio.toFixed(2) + " : " + (100-gcRatio).toFixed(2) + "<br><br>");

  var lastCodon = codons[(codons.length-1)];
  if (stop.length === 0){
    if (lastCodon.length < 3){
    $("#dnaSeqWarn").append("Note: the last codon only has " + lastCodon.length + " base(s). This codon was removed prior to sequence translation.<br><br>");
    translateProt(codons)
    } else{
    translateProt(codons);
      };
};

// Section to deal with internal stop codons
if (stop.length > 0 && stop[0] < codons.length){
    if (stop.length === 1){
      $("#dnaSeqWarn").append("Warning! Your sequence contains a stop codon (bases " + ((stop[0]*3)-2) + " to " + (stop[0]*3) + ") prior to the final codon.")
      $("#dnaSeqWarn").append("<br>All bases after this point were removed prior to translation and analysis below.<br><br>");
    } else if (stop.length > 1){
      $("#dnaSeqWarn").append("Warning! Your sequence contains multiple stop codons before the final codon.")
      $("#dnaSeqWarn").append("<br>All bases after the first of these (bases " + ((stop[0]*3)-2) + " to " + (stop[0]*3) + ") were removed prior to translation and analysis below.<br><br>")
      };
var truncSequence = codons.slice(0,stop[0]-1);
translateProt(truncSequence);
};

if (stop.length < 2 && (lastCodon === "TGA" || lastCodon === "TAA" || lastCodon === "TAG")){
  var codonsLastStop = codons.slice(0,(codons.length-1));
  $("#dnaSeqComp").append("Note: There is a stop codon at the end of your sequence that has been truncated.<br><br>")
  translateProt(codonsLastStop)
};
};

// Function to translate protein
function translateProt(codons){
  var codonTable = {
  AAA: "K",AAC: "N",AAG: "K",AAT: "N",
  ACA: "T",ACC: "T",ACG: "T",ACT: "T",
  AGA: "R",AGC: "S",AGG: "R",AGT: "S",
  ATA: "I",ATC: "I",ATG: "M",ATT: "I",

  CAA: "Q",CAC: "H",CAG: "Q",CAT: "H",
  CCA: "P",CCC: "P",CCG: "P",CCT: "P",
  CGA: "R",CGC: "R",CGG: "R",CGT: "R",
  CTA: "L",CTC: "L",CTG: "L",CTT: "L",

  TAA: ".",TAC: "Y",TAG: ".",TAT: "Y",
  TCA: "S",TCC: "S",TCG: "S",TCT: "S",
  TGA: ".",TGC: "C",TGG: "W",TGT: "C",
  TTA: "L",TTC: "F",TTG: "L",TTT: "F",

  GAA: "E",GAC: "D",GAG: "E",GAT: "D",
  GCA: "A",GCC: "A",GCG: "A",GCT: "A",
  GGA: "G",GGC: "G",GGG: "G",GGT: "G",
  GTA: "V",GTC: "V",GTG: "V",GTT: "V",
  };

var protSeq = [];
for (var i=0; i < codons.length; i++){
  var codon = codons[i];
  for (var cod in codonTable){
     if (cod == codon){
       var res = codonTable[cod];
       protSeq.push(res);
     };
  };
};

var protein = protSeq.join("");
var proteinLength = protein.length;
$("#proteinSeqLength").append("<strong><span>PROTEIN ANALYSIS</span><br>Sequence</strong><br>The translated protein sequence has " + proteinLength + " amino acids:<br>");
aminoComp(protein, proteinLength);
};
