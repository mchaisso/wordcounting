#ifndef BIT_NUC_VECTOR_H_
#define BIT_NUC_VECTOR_H_

#include <limits.h>
#include <string>
#include <vector>

#include "FASTASequence.h"
#include "SG_Types.h"
#include "ByteReverseComplement.h"

const Nuc AdjacencyToBinaryNucleotide[] = {4, 0, 1, 4, 2, 4, 4, 4, 3};
const          Word WordLowerTwo = 0x3;
const unsigned char ByteUpperTwo = 0xc0;
const unsigned char ByteLowerTwo = 0x3;
const Word FullMaskOn = ULONG_MAX;
const Word NucsPerWord = 32;


#define WordMaskOn(pos) ((WordLowerTwo) << (pos))
#define ShiftLeft(nuc, pos) ((nuc) << (pos))
#define ShiftRight(nuc, pos) ((nuc) >> (pos))
#define WordMaskOff(pos)  (~(WordMaskOn(pos)))

using namespace std;

class BitNucVector {
public:
	//
	// Implement a binary encoded nucleotide sequence as a series of words, 16
	//
	
	static Word leftBoundaryMask;
	static Word maxWords;
	static Word numWords;
	static int  k;

	void ToString(string &str, int pk=0) const {
		if (pk == 0) { pk = k; }
		str.resize(pk);
		int i;
		for (i = 0; i < pk; i++) {
			Nuc n;
			Get(i, n);
			str[i] = TwoBitToAscii[n];
		}
	}

	string ToString() const {
		string res;
		ToString(res);
		return res;
	}

	static void Initialize(int _k) {
		k = _k;
		leftBoundaryMask = FullMaskOn;
		int i;
		for (i = 0; i < k % NucsPerWord; i++) {
			leftBoundaryMask <<= 2L;
		}
		leftBoundaryMask = ~leftBoundaryMask;
	}
		
	void Initialize() {
/*
		numWords = 1;
		data = 0;
				for (i = 0; i < numWords; i++) {
			data[i] = 0;
			}*/
	}

	BitNucVector() {
		Initialize();
	}

	BitNucVector(int _k) {
		//
		// This mask is used when clearing upper bits when shifting left.
		//
		k = _k;
		Initialize();
	}
	
	//	Word data[MaxLength/NucsPerWord + ((MaxLength % NucsPerWord) && 1)];
	Word data;
	
	//
	// Set two bits at a particular position
	//
	void Set(int pos, Nuc nuc) {
		int wordIndex = numWords - pos / NucsPerWord - 1;
		Word bitOffset = (pos % NucsPerWord) * 2;
		Word bitMaskOut = WordMaskOff(bitOffset);
		Word wordNuc = nuc;
		Word bitMaskIn = wordNuc << bitOffset;
		//		data[wordIndex] = (data[wordIndex] & bitMaskOut) + bitMaskIn;
		data = (data & bitMaskOut) + bitMaskIn;
	}

	//
	// Get the two bits at a position.
	//
	void Get(int pos, Nuc &nuc) const {
		int wordIndex = numWords - pos / NucsPerWord - 1;
		Word bitOffset = (pos % NucsPerWord) * 2;
		Word bitMask   = WordMaskOn(bitOffset);
		nuc = ((data & bitMask) >> bitOffset) & ByteLowerTwo;		
		//		nuc = ((data[wordIndex] & bitMask) >> bitOffset) & ByteLowerTwo;

	}

	void ShiftLeftTwo() {
		int i;
		data = data << 2L;
		/*
		for (i = 0; i < numWords; i++) {
			data[i] = data[i] << 2;
		}
		data[0] &= leftBoundaryMask;*/
	}

	void ShiftRightTwo() {
		int i;
		Word prevSpill = 0, curSpill;
		//
		// Shift by two, carrying over the out shifted right two bits each
		// new word.
		//
		data = data >> 2L;
		/*
		for (i = 0; i < NumWords(); i++) {
			curSpill  = data[i] & 0x3;
			prevSpill = prevSpill << (Word)((NucsPerWord - 1) * 2);
			data[i]   = data[i] >> 2 + prevSpill;
			prevSpill = curSpill;
			}*/
	}

	void SetReverseComplement(BitNucVector &dest, int kp=0) const {
//#ifdef BYTE_REVERSE_COMPLEMENT
		Word v0, v1, v2, v3, v4, v5, v6, v7;
		if (kp == 0) {
			kp = k;
		}
		v0 = ByteReverseComplement[ data & 0xff];
		v1 = ByteReverseComplement[(data >> 8) & 0xff];
		v2 = ByteReverseComplement[(data >> 16) & 0xff];
		v3 = ByteReverseComplement[(data >> 24) & 0xff];
		v4 = ByteReverseComplement[(data >> 32) & 0xff];
		v5 = ByteReverseComplement[(data >> 40) & 0xff];
		v6 = ByteReverseComplement[(data >> 48) & 0xff];
		v7 = ByteReverseComplement[(data >> 56) & 0xff];
		dest.data = (v0 << 56) 
			+ (v1 << 48)
			+ (v2 << 40)
			+ (v3 << 32) 
			+ (v4 << 24)
			+ (v5 << 16)
			+ (v6 << 8)
			+ (v7);
		dest.data = dest.data >> (64 - 2*kp);
/*#else
		dest.data = 0;
		int i;
		Word complement = ~data;
		assert(k < 32);
		for (i = 0; i < k; i++) {
			dest.data = dest.data << 2L;
			dest.data += (complement & WordLowerTwo);
			complement = complement >> 2L;
		}
#endif
*/
	}

	void Clear() {
		int i;
		data = 0;
		/*		for (i = 0; i < numWords; i++) {
			data[i] = 0;
			}*/
	}
	
	bool InitializeFromString(unsigned char *seq, int length) {
		int i;
		Clear();
		assert(length <= k);
		for (i = 0; i < length; i++) {
			if (TwoBit[seq[i]] > 3)
				return false;
			Set(i, TwoBit[seq[i]]);
		}
		return true;
	}

	bool InitializeFromSequence(DNASequence &seq, int start=0, int length=-1) {
		int i;
		int end = seq.length;
		if (length != -1) {
			end = start + length;
		}
		return InitializeFromString(&seq.seq[start], end - start);
	}

	int operator<(const BitNucVector &rhs) const {
		int i;
		return data < rhs.data;
		/*		for (i = numWords - 1; i >= 0; i--) {
			if (data[i] != rhs.data[i]) {
				return data[i] < rhs.data[i];
			}
			}*/
		return 0;
	}

	int operator==(const BitNucVector &rhs) const {
		int i;
		return data == rhs.data;
		/*
		for (i = numWords - 1; i >= 0; i--) {
			if (data[i] != rhs.data[i]) {
				return false;
			}
			}
		return true;
		*/
	}
		
	BitNucVector &operator=(const BitNucVector &rhs) {
		data = rhs.data;
		//		memcpy(data, rhs.data, sizeof(data));
		return *this;
	}

	void Print(ostream &out) const {
		int i;
		out << (Word) data << endl;
		/*		for (i = 0; i < numWords; i++) {
			out << (Word) data[i] << " ";
		}
		out << endl;*/
	}

	void ShiftOneNucleotide(unsigned char nuc) {
		Word temp = nuc;
		temp = temp << ((k-1)*2);
		data = data >> 2L;
		data += temp;
	}

	void SetKey(BitNucVector &key, int kp=0) const {
		BitNucVector rc;
		SetReverseComplement(rc,kp);

		if (*this < rc) {
			key.data = data;
		}
		else {
			key.data = rc.data;
		}
	}

	Strand GetOrientation() {
		BitNucVector rc;
		SetKey(rc);
		return data == rc.data ? Forward : Reverse;
	}

};

Word data(BitNucVector &v) {
	return v.data;
}


ostream &operator<<(ostream &out, const BitNucVector &b) {
	out << b.data;
	return out;
}

Word BitNucVector::leftBoundaryMask = 0L;
Word BitNucVector::numWords = 0L;

int  BitNucVector::k = 0L;
typedef BitNucVector Tuple;

int size(vector<BitNucVector> &v) {
	return v.size();
}
#endif
