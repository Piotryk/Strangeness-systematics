//dopasowanie funkcji 2d z policzeniem deltaf
#include <fstream>

pair<double,double> DoInterpolateYield ( TF2& fun, int FunType, TMatrixDSym& cov, double A, double Tb, double C, int verbose ) ;


int dop1kpSqSraw4 () {

  ifstream data ("dane1KaonPlus_sqS_tgl.dat");
//ifstream data ("dane1KaonPlus_sqS_tgl_AuAu10.7GeV.dat");
  if (! data.is_open() ) {
    cout << "<E> Input file not found\n";
    return -1;
  }
  ofstream wynik ("wDopasowanie1KPlus.dat");

//Definicje
	//Apart - x, Tbeam- y, multi - 2
	//Apart - alfa, Tbeam - beta
	//N - 0, Apart - 1, Tbeam - 2
	//N - 0, alfa - 1, beta - 2

  Double_t sqS=0., Apart=0., multi=0., errSqS=0., errApart=0., errmulti=0.;
//	Int_t n = liczba_linii;
  Double_t arSqS  [100];
  Double_t arApart[100];
  Double_t armulti[100];
  Double_t arerrSqS  [100];
  Double_t arerrApart[100];
  Double_t arerrmulti[100];
  Int_t i = 0;

/* Wczytanie danych
 
   Tb  | 1.0   | 1.06  | 1.135 | 1.20  | 1.23  | 1.50  | 1.80  | 1.93  |
   sqS | 2.324 | 2.348 | 2.378 | 2.403 | 2.415 | 2.518 | 2.627 | 2.673 |
*/
  while (data >> sqS) 
  {
    data >> Apart >> multi >> errSqS >> errApart >> errmulti;
//    if (! (sqS >= 2.66 && sqS<=2.7 ) ) continue;
//  if (sqS == 2.627 && Apart == 18.11 ) continue;
    arSqS[i] = sqS;
    arApart[i] = Apart;
    armulti[i] = multi;
    arerrSqS[i] = errSqS;
    arerrApart[i] = errApart;
    arerrmulti[i] = errmulti;
    cout << i+1 <<'\t'<< sqS <<'\t'<< Apart << "\t +- " << errApart << "\t P = "
         << multi <<"\t +- " << errmulti << endl;
    i++;
  }

  Int_t liczba_linii = i;


// Deklaracja funkcji fitującej

  int FunType;
  double CONST = 0. ;
/*
  TF2 fit ("fit", "pow(10.,[0]) * x^[1] * y^[2]", 2., 320., 0.5, 3.2 );
  fit.SetParNames ( "N", "a", "b" );
  FunType = 1;
  fit.SetParameter(0, -13. ); // 3e-10);
  fit.SetParameter(1, 1.2 );
  fit.SetParameter(2, 24. );
  fit.SetParLimits(0, -40., -1.);
  fit.SetParLimits(1, 0.9, 2.2);
  fit.SetParLimits(2, 1e-5, 200.);
*/
  CONST = 1080. ;

  TF2 fit ("fit", "[0] * x^[1] * exp( - ([3] * y)^[2] )", 2., 320., 0.5, 3.2);
//  TF2 fit ("fit", "[0] * x^[1] * exp( - [3] * y^[2] )", 2., 320., 0.5, 3.2);

  fit.SetParNames ( "N", "a", "b", "C" );
  FunType = 2;
  fit.SetParameter (0,  3e-3 );
  fit.SetParameter (1,  1.32 );
  fit.SetParameter (2, -6.22 );
  fit.SetParameter (3,  0.323);
  
  fit.SetParLimits (0,  1e-5 , 0.1 );
  fit.SetParLimits (1,  0.6  , 2.0 );
  fit.SetParLimits (2, -9.   , 2.  );
  fit.SetParLimits (3,   0.04 , 0.6 );
/*
  CONST = 1080. ;

  TF2 fit ("fit", Form ("[0] * x^[1] * exp( - %.2f * y^[2] )", CONST) , 2., 320., 0.5, 3.2);
//TF2 fit ("fit", Form ("10^([0]) * x^[1] * exp( - %.2f * y^[2] )", CONST) , 2., 320., 0.5, 3.2);

  fit.SetParNames ( "N", "a", "b" );
  FunType = 3;
  fit.SetParameter (0,  3e-3 );
  fit.SetParameter (1,  1.32 );
  fit.SetParameter (2,  -6.22 );

  fit.SetParLimits (0, 1e-8 , 1e3 );
  fit.SetParLimits (1, 0.6, 2.3 );
  fit.SetParLimits (2, -50. , 2. );
*/
/*
  TF2 fit ("fit", "pow(10.,[0]) * x^([1] - y) * y^[2]", 2., 320., 0.5, 3.2 );
  fit.SetParNames ( "N", "a", "b" );
  FunType = 4;
  fit.SetParameter(0, -10. ); // 3e-10);
  fit.SetParameter(1, 3.9 );
  fit.SetParameter(2, 24. );
  fit.SetParLimits(0, -20., -5.);
  fit.SetParLimits(1, 2., 6.);
  fit.SetParLimits(2, 5., 40.);
*/
/*
  CONST = 350. ;

  TF2 fit ("fit", Form ("[0] * x^([1] - y) * exp( - %.2f / y^[2] )", CONST) , 2., 320., 0.5, 3.2);
  fit.SetParNames ( "N", "a", "b" );
  FunType = 5;
  fit.SetParameter (0,  3e-3 );
  fit.SetParameter (1,  3.9 );
  fit.SetParameter (2,  6.22 );

  fit.SetParLimits (0, 0. , 50. );
  fit.SetParLimits (1, 0.5, 5.5 );
  fit.SetParLimits (2, 2. , 10. );
*/

// Deklaracja grafu i ustawienie opcji

  TGraph2DErrors* graph = new TGraph2DErrors (liczba_linii, arApart, arSqS, armulti, arerrApart, arerrSqS, arerrmulti);

  graph->SetTitle ("Multiplicity with fitted function; Apart; s^{1/2}; Multiplicity");
  graph->SetMarkerSize (0.6);
  graph->SetMarkerStyle (20);
  graph->SetMarkerColor (2);

// Fit

  TVirtualFitter::SetMaxIterations (50000);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy (2);

// For the case of fixed sqS -> f = f(Apart) only, 1-dim fitting.

/*
   TGraphErrors* graphAp = new TGraphErrors 
         (liczba_linii, arApart, armulti, arerrApart, arerrmulti);

   TF1* fitAp = new TF1 ("fitap", "[0] * x^[1]", 1. , 320.);
   fitAp->SetParameters ( 2e-6 , 1.4 ) ;
   fitAp->SetParLimits ( 1 , 0.7 , 2.9 );
  
   gStyle->SetLabelSize (0.055, "XY");
   gStyle->SetNdivisions ( 505, "XY");
   
   TCanvas* c0 = new TCanvas ("c0", "", 640, 480);
   c0->SetLeftMargin (0.15);
   c0->SetRightMargin (0.04);
   c0->SetTopMargin (0.05);
   c0->SetBottomMargin (0.15);
   c0->SetLogy (1);
   
   graphAp->SetMinimum (1e-3);
   graphAp->SetMaximum ( 0.1 );
   graphAp->SetTitle ("");
   graphAp->Draw ("AP");
   graphAp->Fit ( fitAp , "SEM");
   graphAp->Fit ( fitAp , "SEM");
   Double_t alpha = fitAp->GetParameter (1) , dalpha = fitAp->GetParError (1), 
            chiAp = fitAp->GetChisquare ()  , NDFAp = fitAp->GetNDF();
   cout << "NDF = " << NDFAp << endl;
   cout << "Chi2 / NDF = " << chiAp/NDFAp << endl;
   
   TLatex l0;
   l0.SetNDC (1);
   l0.SetTextSize (0.054);
   l0.SetTextFont (42);
   l0.DrawLatex (0.33, 0.22, "T_{Beam} = 1.8A GeV,   #sqrt{s_{NN}} = 2.63 GeV");
   l0.DrawLatex (0.21, 0.85, Form ("#alpha = %4.2f #pm %4.2f", alpha, dalpha) );
   l0.DrawLatex (0.21, 0.76, Form ("#chi^{2}/#nu = %3.1f", chiAp/NDFAp) );
   l0.DrawLatex (0.49, 0.03, "#LT A_{part} #GT_{b}");
   l0.SetTextAngle (90);
   l0.DrawLatex (0.04, 0.49, "Yield");

   return 0;
*/
  cout << "\n* Fitting function (type = " << FunType << ")\n\n";

  graph->Fit (&fit, "SE");
  TFitResultPtr fitRes = graph->Fit (&fit, "SE");
  TMatrixDSym cov = fitRes->GetCovarianceMatrix();

  cout << "\n* Errors via sqrt( cov[i][i] ) : \n";
  string ParNames[] = { "N", "a" , "b" , "C" };
  for ( int p = 0 ; p <= 3 ; p++ ) {
    cout << "d( " << ParNames[p] << " ) = " << sqrt( cov[p][p] ) << endl;
  }

  cout << "\n* COVariance matrix: \n";
  for (int r = 0; r < cov.GetNrows() ; r++) {
    for (int c = 0; c < cov.GetNcols() ; c++)
      cout << setw(16) << cov[r][c] << flush;
    cout << endl;
  }

  cout << "\n* CORrelation matrix: \n";
  for (int r = 0; r < cov.GetNrows() ; r++) {
    for (int c = 0; c < cov.GetNrows() ; c++) 
      cout << setw(16) << cov[r][c] / sqrt(cov[r][r] * cov[c][c]) << flush;
    cout << endl;
  }
//return 0;
// LEK na to, że funkcja fit się nie rysuje ¯\_(ツ)_/¯

  TF2* fit2 = (TF2*)graph->FindObject("fit");
  fit2->SetLineColor(1);
  fit2->SetLineWidth(1);

  Double_t chi = fit.GetChisquare();
  Double_t NDF = fit.GetNDF();
  i = 0;

  cout << endl << "NDF = " << NDF << endl;
  cout << "Chi2 / NDF = " << chi/NDF << "\n\n";

  double Tb, mn = 0.9385, SqrtS, A;
  double A1_23[] = { 156.03 , 111.30 , 79.09 , 55.09 }; // For Ag+Ag
  double A1_58[] = { 159.50 , 114.07 , 80.96 , 56.24 }; // For Ag+Ag
  double A2_91[] = { 311.0  , 161.4  , 57.8  };         // For Au+Au STAR
  
  Tb = 2.1 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A =  9.78;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  Tb = 0.60 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 75.7 ;
  DoInterpolateYield (fit , FunType , cov , A , SqrtS, CONST , true);
  
  cout << endl;
  
  Tb = 1.23 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) ;
  for (int i = 0; i < sizeof( A1_23 ) / sizeof( double ) ; i++ ) 
    DoInterpolateYield (fit , FunType , cov , A1_23[i] , SqrtS, CONST , true);

  cout << endl;
  
  Tb = 1.58 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) ;
  for (int i = 0; i < sizeof( A1_58 ) / sizeof( double ) ; i++ ) 
    DoInterpolateYield (fit , FunType , cov , A1_58[i] , SqrtS, CONST , true);

  cout << endl;

  Tb = 1.23 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 125.8;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  Tb = 1.23 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 149.3;   
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  cout << endl;

  Tb = 1.756,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 37.17;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  cout << endl;
  
  SqrtS = 3.0 ;
  for (int i = 0; i < sizeof( A2_91 ) / sizeof( double ) ; i++ ) 
    DoInterpolateYield (fit , FunType , cov , A2_91[i] , SqrtS, CONST , true);

  cout << endl;
  
  Tb = 10.7 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 341.5;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  A = 100;
  for ( Tb = 2.00 ; Tb <= 11. ; Tb += 2.0 ) {
    cout << "\nT_beam = " << fixed << setprecision (1) << Tb << " AGeV\n" ; 
    SqrtS = sqrt (2.*mn * (2.*mn + Tb)) ;
    DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);
  }

  cout << "\n\n";
  
  Tb = 0.6 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 251.7;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  Tb = 0.6 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 174.8;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

  cout << endl;
  Tb = 0.8 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 290.0;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);
  
  Tb = 0.8 ,  SqrtS = sqrt (2.*mn*(2.*mn + Tb)) , A = 203.6;
  DoInterpolateYield ( fit , FunType , cov , A, SqrtS, CONST , true);

// RYSOWANIE
// Opcje osi grafu - range, ndivisions i położenie labeli

  TAxis* x = fit2->GetXaxis();
  TAxis* y = fit2->GetYaxis();
  TAxis* z = fit2->GetZaxis();
  x->SetNdivisions(505);
  y->SetNdivisions(505);
  z->SetNdivisions(505);

  x->SetTitle("A_{part}");
  x->CenterTitle(true);
  x->SetTitleOffset(1.7);
  y->SetTitle("#sqrt{s_{NN}} [GeV]");
  y->CenterTitle(true);
  y->SetTitleOffset(1.7);
  z->SetTitle("Multiplicity");
  z->CenterTitle(true);
  z->SetTitleOffset(1.7);

//y->SetLimits (0.5, 3.2);
  x->SetLimits (0. , 320);
//Min - max
  graph->SetMinimum ( 0.00001 );
  fit2 ->SetMinimum ( 0.00001 );

  graph->SetMaximum ( 0.8 );
  fit2 ->SetMaximum ( 0.8 );

//Deklaracja TCanvas i ustawienie opcji TCanvas

//Opcje gstyle

  gStyle->SetPalette ( kBird );
  gStyle->SetLabelSize  ( 0.048, "XYZ");
  gStyle->SetNdivisions ( 505 , "XY");

  TCanvas *c1 = new TCanvas("c1", "",   0, 5, 800, 600);
  c1->SetLogz();
  c1->SetLeftMargin   (0.16);
  c1->SetRightMargin  (0.05);
  c1->SetBottomMargin (0.13);
  c1->SetTopMargin    (0.06);

//camera position
  c1->cd ();
  c1->SetTheta (13.8);
  c1->SetPhi   (50.4);

//Rysowanie grafu 3D

  fit2->SetTitle ("");
  fit2->SetNpx (500);
  fit2->SetNpy (500);
//fit2->SetRange (2., 2.2, 10., 2.7);

  fit2->Draw ("surf2");

  graph->Draw ("same err p");

  TLatex l;
  l.SetNDC ( 1 );
  l.SetTextFont ( 42 );
  l.SetTextSize ( 0.078 );
  l.DrawLatex (0.056, 0.89 , "K^{+}" );
  l.SetTextSize ( 0.054 );
  l.SetTextAngle ( 90 );
  l.DrawLatex ( 0.045, 0.465, "Yield" );
  l.SetTextAngle ( -8.2 );
  l.DrawLatex ( 0.25, 0.07, "#sqrt{s}_{NN} [GeV]" );
  l.SetTextAngle ( 12.2 );
  l.DrawLatex ( 0.78, 0.06, "#LTA_{part}#GT_{b}" );

  /*  Drawing residuals scaled with errors, 
      (P_data - P_fit) / sqrt(dP_data^2 - dP_fit^2)     */

//  TH2F* hdiff = new TH2F ("hdiff", "", 80, -0.5, 349.5, 100, 2.1 , 2.7);
  Double_t P, dP, PdevMin = 1e99, PdevMax = 1e-99;
  Double_t Pstdev[100];
  cout << endl;
  
  for (int i = 0; i < liczba_linii; i++) 
  {
    tie (P, dP) = DoInterpolateYield ( fit , FunType , cov , 
                                       arApart[i], arSqS[i], CONST, false);
                                       
    if (arerrmulti[i] < dP) {
      cerr << "<W> Event" << i << " : [dY_data , dY_fit] = [" 
           << arerrmulti[i] << " , " << dP << "]\n" ;
      Pstdev[i] = (armulti[i] - P) / dP ;
    } else {
      Pstdev[i] = (armulti[i] - P) / sqrt (arerrmulti[i]*arerrmulti[i] - dP*dP) ;
    }
    PdevMin = min (PdevMin, Pstdev[i]);
    PdevMax = max (PdevMax, Pstdev[i]);
    cout << "| For Sqrt(s) = " << arSqS[i] << " and Apart = " << arApart[i] 
         << " ,  P_stdev = " << Pstdev[i] << endl;
         
//  hdiff->Fill ( arApart[i] , arSqS[i] , Pstdev ) ;
  }

  cout << " [Min : Max] deviations = [ " << PdevMin << " : " << PdevMax << " ]\n";

  gStyle->SetOptStat (0);
//  gStyle->SetPalette ( 87 , 0 );
  gStyle->SetLabelSize (0.06, "XY");
  gStyle->SetLabelSize (0.052, "Z");
  TCanvas* c2 = new TCanvas ("c2", "", 640, 480);
  c2->SetRightMargin  (0.10);
  c2->SetLeftMargin   (0.17);
  c2->SetTopMargin    (0.05);
  c2->SetBottomMargin (0.18);
  c2->SetTheta (90);
  c2->SetPhi   (0.001);
  c2->SetLogx  (1);
  TGraph2D* grResiduals = new TGraph2D (liczba_linii, arApart, arSqS, Pstdev);
  grResiduals->SetName ("grres");
  grResiduals->SetTitle ("");

  grResiduals->SetMinimum (-9.25);
  grResiduals->SetMaximum ( 3.6);
  grResiduals->SetMarkerSize ( 1.5 );
  grResiduals->SetMarkerStyle (20);
  grResiduals->GetXaxis()->SetLabelOffset( 0.03 ) ;
  grResiduals->GetYaxis()->SetLabelOffset( -0.01 ) ;
  grResiduals->GetYaxis()->SetLimits (2.10, 2.76);
  grResiduals->GetXaxis()->SetLimits (4. , 380);
  grResiduals->Draw ("pcolz");

  TLatex l2;
  l2.SetNDC (1);
  l2.SetTextFont ( 42 );
  l2.SetTextSize ( 0.078 );
  l2.DrawLatex   ( 0.78 , 0.84 , "K^{+}" );
  l2.SetTextSize ( 0.058);
  l2.DrawLatex   ( 0.48, 0.03, "#LTA_{part}#GT_{b}" );
  l2.SetTextAngle( 90. );
  l2.DrawLatex   ( 0.04, 0.43 , "#sqrt{s}_{NN} [GeV]" );

  return 42;
}




pair<double,double> DoInterpolateYield (TF2& fun, int FunType, TMatrixDSym& cov, double A, double SqS, double C = 0. , int verbose = 1)
{
  double P = fun.Eval (A, SqS) , dP = 0.; 
  double Paver, Prms , Ntry, PmeanFit ;

  if (verbose >= 2)
    cout << "+--> Yield from function: " << P << endl;

  ROOT::Math::GSLRandomEngine rnd;
  rnd.Initialize();

  TH1F* hp;
  TH2F* hpN, * hpA , * hpB , * hpC , * hNB , * hNC ;
  
  if (FunType == 1) {
    double N = fun.GetParameter (0) , dN = fun.GetParError (0) , 
           a = fun.GetParameter (1) , da = fun.GetParError (1) ,
           b = fun.GetParameter (2) , db = fun.GetParError (2) ;

    double dPdN = log(10.) * pow(10., N),
           dPda = pow(10., N) * (log(A) * pow(A, a)) * pow(SqS, b),
           dPdb = pow(10., N) * pow(A, a) * (log(SqS) * pow(SqS, b)) ;
                             
    double covNa = cov[0][1] , covNb = cov[0][2] , covab = cov[1][2] ;
  
    dP = pow( dPdN * dN, 2) + pow( dPda * da, 2) + pow( dPdb * db, 2)
            + 2. * dPdN * dPda * covNa + 2. * dPdN * dPdb * covNb + 2. * dPda * dPdb * covab;
    dP = sqrt( dP );
  }
  else if (FunType == 2) 
  {
//  verbose = 2;

    double N = fun.GetParameter (0) , dN = fun.GetParError (0) , 
           a = fun.GetParameter (1) , da = fun.GetParError (1) ,
           b = fun.GetParameter (2) , db = fun.GetParError (2) ,
           C = fun.GetParameter (3) , dC = fun.GetParError (3) ;
/*
N =  2.96e-3 , dN = 9.010e-3 ;
a =  1.316   , da = 2.143e-2 ;
b = -6.184   , db = 0.4891   ;
C =  0.323   , dC = 1.010e-2 ;

cov[0][0] = dN*dN ;
cov[1][1] = da*da ;
cov[2][2] = db*db ;
cov[3][3] = dC*dC ;
cov[0][1] = -5.927e-07, cov[0][2] = 4.177e-4 , cov[0][3] = -8.713e-6 ;
                        cov[1][2] = 2.513e-3 , cov[1][3] = -5.201e-5 ;
                                               cov[2][3] = -4.921e-3 ;
cov[1][0] = cov[0][1], cov[2][0] = cov[0][2], cov[3][0] = cov[0][3] ;
                       cov[2][1] = cov[1][2], cov[3][1] = cov[1][3] ;
                                              cov[3][2] = cov[2][3] ;
*/
                                    
    double dPdN =     pow(A,a)           * exp( - pow(C * SqS, b) ) ,
           dPda = N *(pow(A,a) * log(A)) * exp( - pow(C * SqS, b) ) , 
           dPdb = N * pow(A,a) * exp( -pow(C*SqS, b) ) * (-pow(C*SqS, b) * log(C*SqS) ) ,
           dPdC = N * pow(A,a) * exp( -pow(C*SqS, b) ) * (-b * pow(C*SqS, b-1)) * SqS ;

    double covNa = cov[0][1] , covNb = cov[0][2] , covNC = cov[0][3] ,
                               covab = cov[1][2] , covaC = cov[1][3] ,
                                                   covbC = cov[2][3] ;

    dP = pow( dPdN * dN, 2) + pow( dPda * da, 2) + pow( dPdb * db, 2) + pow( dPdC * dC, 2) 
       + 2. * dPdN * dPda * covNa + 2. * dPdN * dPdb * covNb + 2. * dPdN * dPdC * covNC
                                  + 2. * dPda * dPdb * covab + 2. * dPda * dPdC * covaC
                                                             + 2. * dPdb * dPdC * covbC ;
    dP = sqrt( dP );
    if (verbose >= 2) 
      cout << "\n+--> Analytical method : dP = " << dP << "\n";

    const int dim = 4;
    double*   pars = fun.GetParameters ();
    double genpars[dim] , covMVA[dim * dim] ;

    int indexMVA = 0;
    for (int r = 0; r < dim ; r++)
      for (int c = 0; c < dim ; c++) 
        covMVA[ indexMVA++ ] = cov[r][c] ;
    
    TF2* funCopy = (TF2*) fun.Clone ("funcopy");

    // Estimation of Yield's Average and RMS
    
    Paver = 0., Prms = 0. ;
    for (Ntry = 0; Ntry < 1000 ; Ntry++) {
      rnd.GaussianND ( dim , pars , covMVA , genpars );
      funCopy->SetParameters ( genpars );
      Paver += funCopy->Eval (A, SqS) ; 
      Prms  += pow( funCopy->Eval (A, SqS) , 2. ) ;
    }
    Paver /= Ntry ;
    Prms = sqrt( (Prms - Ntry * Paver*Paver ) / (Ntry - 1) );

    hp  = new TH1F ("hp" , "", 400, Paver- 4.*Prms, Paver+ 4.*Prms);
    hpN = new TH2F ("hpn", "", 100, Paver- 4.*Prms, Paver+ 4.*Prms, 100, 0., pars[0] * 2);
    hpA = new TH2F ("hpa", "", 100, Paver- 4.*Prms, Paver+ 4.*Prms, 100, pars[1] * 0.8, pars[1] * 1.2);
    hpB = new TH2F ("hpb", "", 100, Paver- 4.*Prms, Paver+ 4.*Prms, 100, 0. , pars[2] * 2);
    hpC = new TH2F ("hpc", "", 100, Paver- 4.*Prms, Paver+ 4.*Prms, 100, 0., pars[3] * 2);
    hNB = new TH2F ("hnb", "", 100, -pars[0], pars[0]*3. , 100, -pars[2], pars[2]*3. );
    hNC = new TH2F ("hnc", "", 100, -pars[0], pars[0]*3. , 100, pars[3] * 0.7, pars[3]*1.3 );

    for (Ntry = 0; Ntry < 5e4; Ntry++) {
      rnd.GaussianND ( dim , pars , covMVA , genpars );

      funCopy->SetParameters ( genpars );
      hNB->Fill ( genpars[0] , genpars[2] );
      hNC->Fill ( genpars[0] , genpars[3] );

      hp ->Fill ( funCopy->Eval (A, SqS) );
      hpN->Fill ( funCopy->Eval (A, SqS) , genpars[0] );
      hpA->Fill ( funCopy->Eval (A, SqS) , genpars[1] );
      hpB->Fill ( funCopy->Eval (A, SqS) , genpars[2] );
      hpC->Fill ( funCopy->Eval (A, SqS) , genpars[3] );

      if (verbose >= 3) {
        if (Ntry < 4) {
          for (int i = 0; i < dim; i++)
            cout << setw(12) << genpars[i] ;
            cout << " -> " << funCopy->Eval (A, SqS) << endl;
        }
      }
    }

    double Plt , Prt , hpmax = hp->GetMaximum() ;
    int iLt = 1 , iRt = 1 ;
    double Ytot = hp->Integral( 1, hp->GetNbinsX() ), dPprob ;

    for ( int iScan = 1 ; iScan < hp->GetNbinsX() ; iScan++ ) {
      if ( hp->Integral (1, iScan) < Ytot * 0.5 * (1.0 - 0.683) ) iLt++ ;
      if ( hp->Integral (1, iScan) < Ytot * 0.5 * (1.0 + 0.683) ) iRt++ ;
    }
    Plt = hp->GetBinCenter ( iLt );
    Prt = hp->GetBinCenter ( iRt );
    dPprob = 0.5 * ( Prt - Plt );
    if (iLt == 0 || iRt >= hp->GetNbinsX()) {
      cout << "\n<W-DoInterpolateYield> P-range scan reached limit." ;
      cout << " [iLt : iRt] = [" << iLt << " : " << iRt << "]\n";
      cin.ignore();
    }
    if (verbose >= 2) {
      cout << "\n+--> Method via 68.3 % : dP = " << dPprob 
           << "  within Prob. = " << hp->Integral ( iLt, iRt ) / Ytot << "\n\n";
//exit (0);

      TCanvas* cx = new TCanvas ("cx", "", 640, 480);
      cout << "Showing hNB ... ";
      hNB->Draw("colz");
      cx->Update();
      cin.ignore();
      cout << "Showing hNC ... ";
      hNC->Draw("colz");
      cx->Update();
      cin.ignore();
      cout << "Showing hpN ... ";
      hpN->Draw("colz");
      cx->Update();
      cin.ignore();
      cout << "Showing hpA ... ";
      hpA->Draw("colz");
      cx->Update();
      cin.ignore();
      cout << "Showing hpB ... ";
      hpB->Draw("colz");
      cx->Update();
      cin.ignore();
      cout << "Showing hpC ... ";
      hpC->Draw("colz");
      cx->Update();
      cin.ignore();
      cout << "Showing hp (YIELD) ... ";
      hp->Draw();
      TLine ll (Plt, 0., Plt, hpmax/7. );
      TLine lr (Prt, 0., Prt, hpmax/7. );
      ll.SetLineColor (4) ;
      lr.SetLineColor (4) ;
      ll.Draw ("same");
      lr.Draw ("same");
      cx->Update();
      cin.ignore();
      delete cx;
    }
    delete funCopy;
    delete hNB;
    delete hNC;
    delete hp;
    delete hpN;
    delete hpA;
    delete hpB;
    delete hpC;

    dP = dPprob ;
  }
  else if (FunType == 3) 
  {
    double N = fun.GetParameter (0) , dN = fun.GetParError (0) , 
           a = fun.GetParameter (1) , da = fun.GetParError (1) ,
           b = fun.GetParameter (2) , db = fun.GetParError (2) ;

    double dPdN =     pow(A,a)           * exp( -pow(C* SqS, b) ) ,
           dPda = N *(pow(A,a) * log(A)) * exp( -pow(C* SqS, b) ) , 
           dPdb = N * pow(A,a)           * exp( -pow(C* SqS, b) ) * (-pow(C*SqS, b) * log(C*SqS));

    double covNa = cov[0][1] , covNb = cov[0][2] ,
                               covab = cov[1][2] ;

    dP = pow( dPdN * dN, 2) + pow( dPda * da, 2) + pow( dPdb * db, 2) 
       + 2. * dPdN * dPda * covNa + 2. * dPdN * dPdb * covNb 
                                  + 2. * dPda * dPdb * covab ;
    dP = sqrt( dP );
    
    const int dim = 3;
    double*   pars = fun.GetParameters ();
    double genpars[dim] , covMVA[dim * dim] ;
    PmeanFit = fun.Eval (A, SqS) ;
    Paver = 0., Prms = 0. ;
    
    int indexMVA = 0;
    for (int r = 0; r < dim ; r++)
      for (int c = 0; c < dim ; c++) 
        covMVA[ indexMVA++ ] = cov[r][c] ;

    TF2* funCopy = (TF2*) fun.Clone ("funcopy");

    for (Ntry = 0; Ntry < 2e4; Ntry++) {
      rnd.GaussianND ( dim , pars , covMVA , genpars );
      funCopy->SetParameters ( genpars );

      Paver += funCopy->Eval (A, SqS) ; 
      Prms  += pow( funCopy->Eval (A, SqS) , 2. ) ;
    }
    Paver /= Ntry ;
    Prms = sqrt( (Prms - Ntry * Paver*Paver ) / (Ntry - 1) );

    delete funCopy;
    
    if (verbose) {
      cout << fixed << "For Sqrt(s) = " <<  setprecision(5) << left << SqS << 
              " and Apart = "  <<  setw(5) << setprecision(1) << right << A 
           << " , P = " <<  setprecision(8) << P << " +- " << dP << endl;
      cout << fixed << "For Sqrt(s) = " <<  setprecision(5) << left << SqS << 
              " and Apart = "  <<  setw(5) << setprecision(1) << right << A 
           << " , P = " <<  setprecision(7) << PmeanFit << " +- " << Prms << endl;
    }
    return make_pair ( P , Prms ); 
  }

  if (verbose) {
    cout << fixed << "For Sqrt(s) = " <<  setprecision(5) << left << SqS << 
            " and Apart = "  <<  setw(5) << setprecision(1) << right << A 
         << " , P = " <<  setprecision(7) << P << " +- " << dP << endl;
  }
  return make_pair (P, dP);
}
