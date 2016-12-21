LHE_PropagatorRewgt = [
   "Name:SamplePropagator Alias:<Name> PropScheme:CPS hmass:<HMASS> isGen:1 NoBranch:1 isProp:1
   "Name:CPStoBWPropRewgt Alias:<Name> PropScheme:FixedWidth hmass:<HMASS> Options:DivideP=SamplePropagator isGen:1 isProp:1
]
LHE_Probabilities_MCFM = [
   "Name:SampleHypothesisMCFM Alias:<Name> Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=0.9885296538,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1 NoBranch:1",

   "Name:JJEW_SIG_ghv1_1_MCFM Alias:SampleHypothesisMCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1prime2_1E4_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1_prime2=10000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv2_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv4_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghza1prime2_1E4_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs1_prime2=10000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghza2_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghza4_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_gha2_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_gha4_1_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_SIG_ghv1_1_ghv1prime2_25E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=2500,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv1prime2_25E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,2500 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p251_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.25,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,0.25 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.25,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,0.25 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_25E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=2500,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_25E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,2500 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.25,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,0.25 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.25,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,0.25 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.25,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,0.25 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p25_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.25,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p25i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,0.25 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_SIG_ghv1_1_ghv1prime2_50E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=5000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv1prime2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,5000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv1prime2_50E2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=5000,5000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.5,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.5,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_50E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=5000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,5000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_50E2_50E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=5000,5000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.5,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.5,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.5,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p5_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.5,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p5_0p5i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_SIG_ghv1_1_ghv1prime2_75E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=7500,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv1prime2_75E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,7500 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p751_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0.75,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv2_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,0.75 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0.75,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghv4_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,0.75 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_75E2_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=7500,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza1prime2_75E2i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,7500 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0.75,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza2_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,0.75 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0.75,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_ghza4_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,0.75 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0.75,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha2_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,0.75 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p75_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0.75,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_SIG_ghv1_1_gha4_0p75i_MCFM Process:HSMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,0.75 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_BSI_ghv1_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1_prime2=10000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghza1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs1_prime2=10000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghza2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghza4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_gha2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_gha4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_BSI_ghv1_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1prime2_50E2_50E2i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1_prime2=5000,5000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv2_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz2=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv4_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz4=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghza1prime2_50E2_50E2i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs1_prime2=5000,5000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghza2_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs2=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghza4_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghzgs4=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_gha2_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs2=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_gha4_0p5_0p5i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghgsgs4=0.5,0.5 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_BSI_ghv1_1_ghv1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=10000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv1prime2_1E4i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz1_prime2=0,10000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv2_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz2=0,1 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghv4_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghz4=0,1 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza1prime2_1E4_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=10000,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza1prime2_1E4i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs1_prime2=0,10000 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza2_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs2=0,1 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_ghza4_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghzgs4=0,1 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_gha2_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_gha2_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs2=0,1 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_gha4_1_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=1,0 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJEW_BSI_ghv1_1_gha4_i_MCFM Process:bkgZZ_SMHiggs Production:JJEW MatrixElement:MCFM Couplings:ghz1=1,0;ghgsgs4=0,1 hmass:125 hwidth:0.00407 Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",

   "Name:JJEW_BKG_MCFM Process:bkgZZ Production:JJEW MatrixElement:MCFM Options:DivideP=SampleHypothesisMCFM Cluster:BestLOAssociatedVBF isGen:1",
   "Name:JJQCD_BKG_MCFM Process:bkgZZ Production:JJQCD MatrixElement:MCFM Options:DivideP=SampleHypothesisMCFM Cluster:CommonLast isGen:1",
]

# Construct the final list
theLHEProbabilities = []
theLHEProbabilities.extend(LHE_PropagatorRewgt)
theLHEProbabilities.extend(LHE_Probabilities_MCFM)

# Append final list
for name in (
             "ZZTree",
             "CRZLLTree",
             "ZZTreelooseEle",
             "CRZLLTreelooseEle",
             "CRZLLTreeZ1RSE",
             "ZZTreetle",
             "CRZLLTreetle",
            ):
    if hasattr(process, name):
        getattr(process, name).lheProbabilities.extend(theLHEProbabilities)
