

all: findUnique guaranteeUnique queryBam queryFasta countRepeatComposition filterlc nloc

BLASR_COMMON=blasr/common
SAMTOOLS=samtools
FLAGS=-g -std=c++98

findUnique: FindUnique.cpp
	g++ $(FLAGS) FindUnique.cpp -I $(BLASR_COMMON) -o findUnique

countRepeatComposition: CountRepeatComposition.cpp
	g++ $(FLAGS) CountRepeatComposition.cpp -I $(BLASR_COMMON) -o countRepeatComposition -static

guaranteeUnique: GuaranteeUnique.cpp
	g++ $(FLAGS) GuaranteeUnique.cpp -I $(BLASR_COMMON) -o guaranteeUnique -lpthread


queryBam: QueryBam.cpp
	g++ $(FLAGS) QueryBam.cpp -I $(PWD) -I $(BLASR_COMMON) -I htslib/include/htslib -L htslib/lib -lhts -lbam -lm -lz -lpthread -o $@ 

queryFasta: QueryFasta.cpp
	g++ $(FLAGS)  QueryFasta.cpp -I $(BLASR_COMMON) -lm -L$(HOMEDIR)/software/lib -lz -lpthread -o $@ -static

filterlc: FilterLowComplexity.cpp
	g++ $(FLAGS) $^ -I $(BLASR_COMMON)  -lm -lz -lpthread -o $@ -static

nloc: PrintNLocations.cpp
	g++ $(FLAGS) $^ -I $(BLASR_COMMON) -lpthread -o $@ -static

