#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

int main () {
	string filename = "";
	string nieistotne = "hehe";
	int ID = 0;
	char first = 'x';
	float parameter = 0.; //parametr zdzeneia
	float totalevent = 0, totalkaonzero = 0, totalkaonplus = 0, totalkaonminus = 0, totalPhiKminus = 0, totalPhiReconstructed = 0, totallambda = 0, totalsigma = 0, totalPiPlus = 0, totalPiMinus = 0, totalPiZero = 0;
	int totalfiles = 0, failed_files = 0, file_fail_flag = 0;
	int tempkaonplus = 0, tempkaonminus = 0, tempkaonzero = 0, tempPhiKminus = 0, tempPhiReconstructed = 0, templambda = 0, tempsigma = 0,  tempPiPlus = 0, tempPiMinus = 0, tempPiZero = 0;
	int dlugosc = 7; //liczba klas centralnosci
	int   nevent				[ dlugosc ]	= {0,0,0,0,0,0,0};
	float kaonzero				[ dlugosc ] = {0,0,0,0,0,0,0};
	float kaonplus				[ dlugosc ] = {0,0,0,0,0,0,0};
	float kaonminus			[ dlugosc ] = {0,0,0,0,0,0,0};
	float PhiKminus			[ dlugosc ]	= {0,0,0,0,0,0,0};
	float PhiAlive				[ dlugosc ]	= {0,0,0,0,0,0,0};
	float PhiReconstructed	[ dlugosc ] = {0,0,0,0,0,0,0};
	float lambda				[ dlugosc ]	= {0,0,0,0,0,0,0};
	float sigma					[ dlugosc ]	= {0,0,0,0,0,0,0};
	float PiPlus				[ dlugosc ]	= {0,0,0,0,0,0,0};
	float PiMinus				[ dlugosc ]	= {0,0,0,0,0,0,0};
	float PiZero				[ dlugosc ]	= {0,0,0,0,0,0,0};
	float podzial				[ dlugosc ]	= {3.8, 5.4, 6.6, 7.6, 8.5, 9.3, 10.1};
	//float podzial				[ dlugosc ]	= {5.0, 5.1, 6.6, 7.6, 8.5, 9.3, 10.101}; // Dla Ar + KCl Wtedy istotne są tylko 2 pierwsze klasy. b = 5.0 -> cent = 34% ; b = 5.1 -> cent = 36%
	int centr			  [ dlugosc + 1 ]	= {0,  10, 	20,  30,  40,  50,  60,  100};
	//int centr			  [ dlugosc + 1 ]	= {0,  34, 	36,  99,  99,  99,  99,  99}; // // Dla Ar + KCl
	int flagover100 = 1, over100 = 0; //zliczanie czastek poza bmax
	int number_of_events_in_file = 0; //liczba eventow w pliku
	int nparticles = 0;
	int parent1 = 0, parent2 = 0;

	ifstream files ("in_files_SMASH.txt");
	ofstream wynik ("krotnosci_smash_rec.dat");
	ofstream rec 	("reconstruction_debug.dat");


	//rekonstrukcja Phi
	bool flagPhiP = 0, flagPhiM = 0;
	float px  = 0, py  = 0, pz  = 0, E  = 0, mass = 0; //assuming GeV everywhere
	float px1 = 0, py1 = 0, pz1 = 0, E1 = 0;	//K+
	float px2 = 0, py2 = 0, pz2 = 0, E2 = 0;	//K-

	while (files >> filename){

		ifstream data (filename.c_str());
		//Pomiń  nagłówek pliku
		for (int asd = 1; asd<=3; asd++) {getline(data, nieistotne);}
	  totalfiles++;

		while (data >> first) {
			number_of_events_in_file++; //liczba eventow w pliku
			//Wczytaj linijkę z numerem eventu
			data >> nieistotne >> nieistotne >> nieistotne >> nparticles;

			for (int i = 1; i <= nparticles; i++) {
				data 	>> nieistotne >> nieistotne >> nieistotne >> nieistotne >> nieistotne
						>> E >> px >> py >> pz
				 		>>ID;
				for (int j = 0; j < 8; j++) {
					data >> nieistotne;
				}
				data >> parent1 >> parent2;

				//ilosc cziastek w jednym evencie
				if ( ID == 321 ) {
					tempkaonplus++;
					if (parent1 == 333 && parent2 == 0) {
						px1 = px;
						py1 = py;
						pz1 = pz;
						E1 = E;
						if (flagPhiP ) {cout<<"ERROR K+:"<< tempkaonplus <<endl;}
						flagPhiP = 1;
					} //Pamiętać o podzieleniu przez 48%
				}
				else if ( ID == -321 ) {
					tempkaonminus++;
					if (parent1 == 333 && parent2 == 0) {
						tempPhiKminus++;
						px2 = px;
						py2 = py;
						pz2 = pz;
						E2 = E;
						if (flagPhiM ) {cout<<"ERROR 2K-"<<endl;}
					flagPhiM = 1;
					}
				}
				else if ( ID == 333  ) {cout<<"Error: Phi przeżyło"<<endl;}
				else if ( ID == 3122 ) { templambda++;		}
				else if ( ID == 3212 ) { tempsigma++;		}
				else if ( ID == 311  ) { tempkaonzero++;	}
				else if ( ID == 211  ) { tempPiPlus++;		}
				else if ( ID == -211 ) { tempPiMinus++;	}
				else if ( ID == 111  ) { tempPiZero++;		}

			}
			//wczytaj ostatnią linijkę eventu wraz z parametretrem zderzenia
			for (int j = 0; j < 6; j++) {
				data >> nieistotne;
			}
			data >> parameter >> nieistotne >> nieistotne;


			//rekonstrukcja phi
			if ( flagPhiP && flagPhiM ) {
				px = px1 + px2;
				py = py1 + py2;
				pz = pz1 + pz2;
				E = E1 + E2;
				mass = pow(E,2) - ( pow(px,2) + pow(py,2) + pow(pz,2) );
				mass = pow(mass, 0.5);
				rec << setw(6) << setprecision(6) << mass << endl;	//debug - zrzuca wszystkie masy z pary K+K-
				if ( mass > 1.005 && mass < 1.045 ) {tempPhiReconstructed++;} //SPRAWDZANIE MASY PHI
			}
			flagPhiP = flagPhiM = 0;


			for (int i = 0; i < dlugosc; i++){
				if (parameter < podzial[i]){
					flagover100 			=  0;
					nevent[i] 				+= 1;
					kaonplus[i] 			+= tempkaonplus;
					kaonminus[i] 			+= tempkaonminus;
					kaonzero[i] 			+= tempkaonzero;
					PhiKminus[i]	 		+= tempPhiKminus;
					PhiReconstructed[i] 	+= tempPhiReconstructed;
					sigma[i] 				+= tempsigma;
					lambda[i] 				+= templambda;
					PiPlus[i] 				+= tempPiPlus;
					PiMinus[i]	 			+= tempPiMinus;
					PiZero[i]				+= tempPiZero;
					break;
				}
			}

			if (flagover100) {over100++;}
			flagover100 = 1;

			//total number
			if(!file_fail_flag) {
				totalkaonplus 			 += tempkaonplus;
				totalkaonminus	 		 += tempkaonminus;
				totalkaonzero 			 += tempkaonzero;
				totalPhiKminus			 += tempPhiKminus;
				totalPhiReconstructed += tempPhiReconstructed;
				totallambda 		 	 += templambda;
				totalsigma	 			 += tempsigma;
				totalPiPlus 			 += tempPiPlus,
				totalPiMinus	 		 += tempPiMinus,
				totalPiZero 			 += tempPiZero;
			}

			//wyczyść temp ID
			tempkaonplus = tempkaonminus = tempkaonzero = tempPhiKminus = tempPhiReconstructed = templambda = tempsigma = tempPiPlus = tempPiMinus = tempPiZero = 0;
		}
		if(!file_fail_flag) {
			cout << "Plik " << filename << " zliczony. eventów było " << number_of_events_in_file << endl;
		}
		totalevent += number_of_events_in_file;
		number_of_events_in_file = 0; //liczba eventow w pliku
		data.close();
	}

	wynik << "liczba plikow = " << totalfiles << endl;
	wynik << "Nevent = " << totalevent << endl;
	wynik << "over100 = " << over100 << endl << endl;
	totalevent = 0;

	wynik << "Liczby zliczeń " << endl;
	wynik << "klasa " <<'\t'<< "Nevents" <<'\t'
			<< "K+" <<'\t'<< "K-" <<'\t'
			<< "PhiKminus" <<'\t'<< "PhiReconstructed" <<'\t'
			<< "Lambda" <<'\t'<< "Sigma0" <<'\t'
			<< "K0" <<'\t'
			<< "PiPlus" <<'\t'<< "PiMinus" <<'\t'<< "PiZero"<<endl;

	for (int i = 0; i < dlugosc; i++){
		wynik << "cent " << centr[i] <<"-"<< centr[i+1] <<'\t'<< nevent[i] <<'\t'
				<< kaonplus[i] <<'\t'
				<< kaonminus[i] <<'\t'
				<< PhiKminus[i] <<'\t'
				<< PhiReconstructed[i] <<'\t'
				<< lambda[i] <<'\t'
				<< sigma[i] <<'\t'
				<< kaonzero[i] <<'\t'
				<< PiPlus[i] <<'\t'
				<< PiMinus[i] <<'\t'
				<< PiZero[i]
				<< endl;
		totalevent += nevent[i];
	}
	wynik << endl << "cent 0-100" <<'\t'<< totalevent <<'\t'
			<< totalkaonplus <<'\t'<< totalkaonminus <<'\t'
			<< totalPhiKminus <<'\t'<< totalPhiReconstructed <<'\t'
			<< totallambda <<'\t'<< totalsigma <<'\t'
			<< totalkaonzero <<'\t'
			<< totalPiPlus <<'\t'<< totalPiMinus <<'\t'<< totalPiZero
			<< endl;

	wynik << endl << endl << "Krotnosci" << endl;
	wynik << "klasa " <<'\t'<< "Nevents" <<'\t'
			<< "K+" <<'\t'<< "K-" <<'\t'
			<< "PhiKminus" <<'\t'<< "PhiReconstructed" <<'\t'
			<< "Lambda" <<'\t'<< "Sigma0" <<'\t'
			<<  "K0" <<'\t'
			<< "PiPlus" <<'\t'<< "PiMinus" <<'\t'<< "PiZero"
			<<endl;

	for (int i = 0; i < dlugosc; i++){
		wynik << setprecision(9) << "cent "<< centr[i] <<"-"<< centr[i+1] <<'\t' << nevent[i] <<'\t'
				<< setprecision(3) << kaonplus[i] / nevent[i] <<'\t'
				<< kaonminus[i] / nevent[i] <<'\t'
				<< PhiKminus[i] / nevent[i] <<'\t'
				<< PhiReconstructed[i] / nevent[i] <<'\t'
				<< lambda[i] / nevent[i] <<'\t'
				<< sigma[i] / nevent[i] <<'\t'
				<< kaonzero[i] / nevent[i] <<'\t'
				<< PiPlus[i] / nevent[i] <<'\t'
				<< PiMinus[i] / nevent[i] <<'\t'
				<< PiZero[i] / nevent[i]<< endl;
	}

	wynik << setprecision(9) << endl << "cent 0-100" <<'\t'<< totalevent <<'\t'
			<< setprecision(3) << totalkaonplus / totalevent <<'\t'
			<< totalkaonminus / totalevent <<'\t'
			<< totalPhiKminus / totalevent <<'\t'
			<< totalPhiReconstructed / totalevent <<'\t'
			<< totallambda / totalevent <<'\t'
			<< totalsigma / totalevent <<'\t'
			<< totalkaonzero / totalevent <<'\t'
			<< totalPiPlus / totalevent <<'\t'
			<< totalPiMinus / totalevent <<'\t'
			<< totalPiZero / totalevent << endl;


	wynik << endl << endl << endl << "phi:" << endl
		   << '\t' << "PhiKminus" <<'\t'<< "PhiReconstructed" << endl;
	for (int i = 0; i < dlugosc; i++){
		wynik << "cent " << centr[i] <<"-"<< centr[i+1] <<'\t'
				<< PhiKminus[i] <<'\t'
				<< PhiReconstructed[i] <<'\t'
				<< endl;
	}


	files.close();
	wynik.close();
	rec.close();
	return 42;
}
