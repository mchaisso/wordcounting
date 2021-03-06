#include "BitNucVector.h"

#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <numeric>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
                         //0, 1, 2, 3, 4, 5, 6, 7, 8
const int BamToTwoBit[] = {0, 0, 1, 0, 2, 0, 0, 0, 3};
const char BamToAscii[] = {0,'A','C',0,'G',0,0,0,'T',0,0,0,0,0,0,0,0,0,0,0,'N'};

sem_t *semreader;
sem_t *semcount;
sem_t *semfastqwriter;

void rand_str(char *dest, size_t length) {
    char charset[] = "0123456789"
                     "abcdefghijklmnopqrstuvwxyz"
                     "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		
    while (length-- > 0) {
        size_t index = (double) rand() / RAND_MAX * (sizeof charset - 1);
        *dest++ = charset[index];
    }
    *dest = '\0';
}

int AdvanceNucleotide(string &seq, int pos, int length, BitNucVector &tuple) {

	if (pos +tuple.k >= length) {
		return pos + tuple.k;
	}
	if (seq[pos] != 'N') {
		tuple.ShiftOneNucleotide(TwoBit[seq[pos+tuple.k]]);
		// 
		// this nucleotide is fine, next iteration will not reset tuple.
    return 0;
	}
	else {
		// signal that must start search past 'N'
		return pos + tuple.k + 1;
	}
}

class Counts {
public:
	ifstream *in, *pairIn;
	ofstream *readsOut;
	bool isFastq;
	bool noDoubleCount;
	map<BitNucVector, vector<int> > *tupleToIndex;
	vector<int> *queryCount;
	int *readIndex ;
	int maxNReads;
	bool interleaved;
	int threadIndex;
	Counts() {
		pairIn = NULL;
		readsOut = NULL;
		interleaved = false;
	}
};

int readNumber = 0;
int nHits = 0;

void CountWords(void *data) {
	map<BitNucVector, vector<int> > &tupleToIndex =  *((Counts*)data)->tupleToIndex;
	vector<int> &queryCount = *((Counts*)data)->queryCount;
	int &readNumber = *((Counts*)data)->readIndex;
	ifstream *in = ((Counts*)data)->in;
	bool isFastq = ((Counts*)data)->isFastq;
	int numSeq = 1;
	int maxNReads = ((Counts*)data)->maxNReads;
	
	ifstream *pairedIn = ((Counts*)data)->pairIn;

	if (((Counts*)data)->interleaved) {
		pairedIn = in;
		numSeq = 2;
	}
	else {
		pairedIn = NULL;
	}

	ofstream *readsOut = ((Counts*)data)->readsOut;
	if (pairedIn != NULL) {
		numSeq = 2;
	}
	bool noDoubleCount = ((Counts*)data)->noDoubleCount;

	map<BitNucVector, vector<int> >::iterator searchResult;
	while (true) {
		int retval;
		//
		// Make sure that only one thread at a time enters the reading section
		//

		retval = sem_wait(semreader);		
		string title, seq, qualtitle, qual;
		string mptitle, mpseq, mpqualtitle, mpqual;

		int totalSize = 0;
		int totalMatches = 0;
		vector<string> seqs;
		vector<string> quals;
		vector<string> titles;

		while (totalSize < 50000000 and (*in)) {
			seq = "";
			if (isFastq) {

				getline(*in, title);
				(*in) >> seq;
				(*in).get();
				getline(*in, qualtitle);
				(*in) >> qual;
				(*in).get();
				seqs.push_back(seq);
				quals.push_back(qual);
				titles.push_back(title);
				totalSize += seq.size();

				if (numSeq == 2) {
					
					getline(*pairedIn, mptitle);
					(*pairedIn) >> mpseq;
					(*pairedIn).get();
					getline(*pairedIn, mpqualtitle);
					(*pairedIn) >> mpqual;
					(*pairedIn).get();
					seqs.push_back(mpseq);
					quals.push_back(qual);
					titles.push_back(mptitle);
					totalSize += mpseq.size();
				}
			}
			else {
				if (numSeq == 2) {
					cout << "ERROR, paired reads only work with fastq for now!" << endl;
					exit(1);
				}
				getline(*in, title);
				while (in->peek() != '>' and in) {

					char p = in->peek();
					//				cout << "got peek: "<< p << endl;
					string line;
					(*in) >> line;
					if (line == "") {
						break;
					}
					assert(line[0] != '>');
					seq += line;
					//				cout << "seq tmp: " << seq << endl;
					if (in->peek() == '\n') {
						in->get();
					}
				}
				seqs.push_back(seq);
				totalSize += seq.size();
			}
			++readNumber;
			if ((readNumber % 1000000 == 0) or (isFastq == false and readNumber  % 10000 == 0)) {
				cerr << readNumber << "\t" << totalSize << "\t" << nHits << endl;
			}


		}
		cerr << "Buffered reading " << totalSize << "\t" << readNumber << endl;
		
		//
		// All done reading, release the thread lock.
		//
		sem_post(semreader);

		
		int seqi;
		bool stop = false;
		bool foundMatch = false;
		for (seqi = 0; seqi < seqs.size() and stop == false; ++seqi) {
			if (numSeq == 2 and seqi % 2 == 0) {
				foundMatch = false;
			}
			if (numSeq == 1) {
				foundMatch = false;
			}
			int seqPos = 0;
			int res;
			BitNucVector tuple;
			bool initializeTuple = true;
			string seq;
			seq = seqs[seqi];
			if (seq.size() < tuple.k) {
				continue;
			}
			do {
				//
				// Store result.
				//
				if (initializeTuple) {
					while (seqPos < seq.size() - tuple.k) {
						if (tuple.InitializeFromString((unsigned char*) &seq[seqPos], tuple.k)) {
							break;
						}
						else {
							seqPos++;
						}
					}
					initializeTuple = false;
				}

				if (seqPos >= seq.size() - tuple.k) {
					break;
				}

				searchResult = tupleToIndex.find(tuple);
				if (searchResult != tupleToIndex.end()) {
					//					sem_wait(semcount);
					int i;
					for (i = 0; i < (*searchResult).second.size(); i++) {
						int index = (*searchResult).second[i];
						queryCount[index] += 1;				
						++nHits;
						foundMatch = true;
						if (noDoubleCount) {
							seqPos = seq.size();
							stop = true;
						}
					}
					//					sem_post(semcount);
				}
				
				// 
				// Move to the next nucleotide;
				//
				int res = AdvanceNucleotide(seq, seqPos, seq.size(), tuple);
				if (res == 0) {
					seqPos++;
				}
				else {
					initializeTuple = true;
				}
			} while (seqPos < seq.size() - tuple.k);
				
			
			if (foundMatch and readsOut != NULL) {
				if (numSeq == 2 and seqi % 2 == 1) {
					sem_wait(semfastqwriter);
					int i;
					for (i = 0; i < numSeq; i++) {
						*readsOut << titles[seqi-1] << endl;
						*readsOut << seqs[seqi-1] << endl;
						*readsOut << qualtitle << endl;
						*readsOut << quals[seqi-1] << endl;
					}
					sem_post(semfastqwriter);			
				}
			}
		}

		if ((*in).good() == false or 
				(isFastq == true and qual == "") or
				(isFastq == false and (title == "" or seq == "")) or
				(maxNReads > 0 and readNumber >= maxNReads)) {
			//			cout << "qual: " << quals[0] << endl << " title " << titles[0] << " qualtitle: " << qualtitle << " seq " << seqs[0] << endl;
			cout << "Finished at read index " << readNumber << endl;
			//		sem_post(semreader);
			return;
		}

	}
}

int main(int argc, char* argv[]) {
	string queryFileName, fastqFileName, outputFileName;
	
	if (argc < 4) {
		cout << endl;
		cout << "Usage: queryFasta file.queries file.fastq outputFile [-single] [-nproc n] [-fasta ] " << endl
				 << "        [-pair paired reads file] [-matches output.fastq]"  << endl; 
		cout << "       file.queries is in the format name seq [pos], and may be constructed " << endl
				 << "       using filterUnique.  Each query sequence must be less than 32 nt long." << endl;
		cout << " Options:" << endl;
		cout << "  -single    Only count the first occurrence of a word (false)." << endl;
		cout << "  -nproc n   Use n threads (1)." << endl;
		cout << "  -fasta     Force FASTA format reading (false)." <<endl;
		cout << "  -pair file Read pair from this file (unpaired)." << endl;
		cout << "  -interleved  Paired sequences are interleaved (no pairing)." << endl;
		cout << "  -matches   Print sequences (and mate-pair if -pair is used) matching " << endl
				 << "             the query table." << endl << endl;
		exit(0);
	}
	queryFileName  = argv[1];
	fastqFileName  = argv[2];
	outputFileName = argv[3];
	int nproc = 6;
	bool forceFasta = false;
	bool noDoubleCount = false;
	string queryMPFileName = "";
	string matchedReadsFileName = "";
	int maxNReads = 0;
	bool interleaved = false;
	if (argc > 4) {
		int argi = 4;
		while (argi < argc) {
			if (strcmp(argv[argi], "-nproc") == 0) {
				nproc = atoi(argv[++argi]);
			}
			else if (strcmp(argv[argi], "-fasta") == 0) {
				forceFasta = true;
			}
			else if (strcmp(argv[argi], "-single") == 0) {
				noDoubleCount = true;
			}
			else if (strcmp(argv[argi], "-pair") == 0) {
				queryMPFileName = argv[++argi];
			}
			else if (strcmp(argv[argi], "-interleaved") == 0) {
				interleaved = true;
			}
			else if (strcmp(argv[argi], "-matches") == 0) {
				matchedReadsFileName = argv[++argi];
			}
			else if (strcmp(argv[argi], "-nreads") == 0) {
				maxNReads = atoi(argv[++argi]);
			}
			++argi;
		}
	}

	ofstream outFile(outputFileName.c_str());
	ifstream queryFile(queryFileName.c_str());
	ifstream fastqFile(fastqFileName.c_str());
	
	ifstream queryMPFile;

	if (queryMPFileName != "") {
		queryMPFile.open(queryMPFileName.c_str());
		interleaved = false;
	}

	ofstream readsOut;
	bool printMatchingReads = false;
	if (matchedReadsFileName != "") {
		readsOut.open(matchedReadsFileName.c_str());
		printMatchingReads = true;
	}

	bool isFastq;
	if (forceFasta) {
		isFastq = false;
	}
	else {
		if (fastqFile.peek() == '>') {
			isFastq = false;
		}
		else if (fastqFile.peek() == '@') {
			isFastq = true;
		}
		else {
			if (fastqFileName.find(".fasta") != fastqFileName.npos) {
				isFastq = false;
			}
			else {
				isFastq = true;
			}
		}
	}
	//
	// Make sure there are no duplicate keys.
	//
	if (not queryFile) {
		cout << "Could not open " << queryFileName << endl; exit(1);
	}
	map<BitNucVector, int> queryList;

	
	vector<int> queryCount;
	vector<string> queryNames;
	vector<int> queryPositions;
	vector<BitNucVector> queryTuples;
	map<BitNucVector, vector<int> > tupleToIndex;
	bool kIsInitialized = false;
	int index = 0;

	while (queryFile) {
		string seq, tupleStr;
		int pos;
		if (! (queryFile >> seq >> tupleStr >> pos)  ) { break; }


		BitNucVector tuple, tuplerc;
		if (kIsInitialized == false) {
			tuple.k = tupleStr.size();
			kIsInitialized = true;
		}
		tuple.InitializeFromString((unsigned char*) tupleStr.c_str(), tuple.k);
		queryCount.push_back(0);
		queryNames.push_back(seq);
		queryTuples.push_back(tuple);
		queryPositions.push_back(pos);
		tupleToIndex[tuple].push_back(index);
		tuple.SetReverseComplement(tuplerc);
		tupleToIndex[tuplerc].push_back(index);
		++index;
		
	}

	map<BitNucVector, int>::iterator searchResult;
	int readNumber = 0;
	
	vector<Counts> counts(nproc);
	int i;
	for (i = 0; i < nproc; i++) {
		counts[i].in = &fastqFile;
		counts[i].interleaved = interleaved;
		counts[i].pairIn = &queryMPFile;
		counts[i].isFastq = isFastq;
		counts[i].tupleToIndex = &tupleToIndex;
		counts[i].queryCount = &queryCount;
		counts[i].readIndex = &readNumber;
		counts[i].noDoubleCount = noDoubleCount;
		if (printMatchingReads) {
			counts[i].readsOut = &readsOut;
		}
		counts[i].maxNReads = maxNReads;
		counts[i].threadIndex = i;
	}
	const int idLen=10;
	char id[idLen+1];
	id[idLen] = '\0';
	srand (time(NULL));
	rand_str(id, idLen);

	string readerName = string("/semreader_") + string(id);
	string countName  = string("/semcount_") + string(id);
	string semfastqwriterName = string("/semfastqwriter_") + string(id);

	semreader     = sem_open(readerName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << readerName << endl;
		exit(1);
	}
	semcount      = sem_open(countName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << countName << endl;
		exit(1);
	}
	semfastqwriter = sem_open(semfastqwriterName.c_str(), O_CREAT, 0644, 1);
	if (semfastqwriter == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << semfastqwriterName << endl;
		exit(1);
	}

	sem_init(semreader, 1, 1);
	sem_init(semcount, 1, 1);
	sem_init(semcount, 1, 1);

	pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
	int t;	

	for (t = 0; t < nproc; t++ ){
		pthread_attr_init(&threadAttr[t]);
	}
	pthread_t *threads = new pthread_t[nproc];

	for (t = 0; t < nproc; t++) {
		pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords, &counts[t]);
	}

	for (t = 0; t < nproc; t++) {
		pthread_join(threads[t], NULL);
	}


	for (i = 0; i < queryNames.size(); i++) {
		string tupleString = queryTuples[i].ToString();
		outFile << queryNames[i] << "\t" << tupleString << "\t" << queryCount[i] << "\t" << queryPositions[i] << endl;
	}

	if (matchedReadsFileName != "") {
		readsOut.close();
	}

	outFile.close();
}


