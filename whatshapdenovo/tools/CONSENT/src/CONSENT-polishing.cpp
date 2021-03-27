#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include "CONSENT-polishing.h"
#include "../BMEAN/bmean.h"
#include "../BMEAN/utils.h"
#include "DBG.h"
#include "correctionAlignment.h"
#include "correctionDBG.h"
#include "correctionMSA.h"
#include "alignmentPiles.h"
#include "alignmentWindows.h"
#include "../CTPL/ctpl_stl.h"

std::mutex outMtx;
robin_hood::unordered_map<std::string, std::vector<bool>> readIndex;
bool doTrimRead = false;

std::pair<std::string, std::string> processContig(std::vector<Overlap>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors,unsigned solidThresh, unsigned windowOverlap, unsigned maxMSA, std::string path, unsigned nbThreads) {
	std::string readId = alignments.begin()->qName;
	robin_hood::unordered_map<std::string, std::string> sequences = getSequencesMap(alignments, readIndex);
	std::vector<std::pair<unsigned, unsigned>> pilesPos = getAlignmentWindowsPositions(alignments.begin()->qLength, alignments, minSupport, maxSupport, windowSize, windowOverlap);
	if (pilesPos.size() == 0) {
		return std::make_pair(readId, "");
	}

	// Compute consensuses for all the piles
	std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>> resCons;
	std::vector<std::string> consensuses(pilesPos.size());
	std::vector<robin_hood::unordered_map<kmer, unsigned>> merCounts(pilesPos.size()); 
	std::vector<std::string> curPile;
	std::vector<std::string> templates(pilesPos.size());

	int poolSize = 100000;
	ctpl::thread_pool myPool(nbThreads);
	int jobsToProcess = pilesPos.size();
	int jobsLoaded = 0;
	int jobsCompleted = 0;

	std::string curTpl;

	// Load the first jobs
	vector<std::future<std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>>>> results(poolSize);
    while (jobsLoaded < poolSize && jobsLoaded < jobsToProcess) {
    	curPile = getAlignmentWindowsSequences(alignments, minSupport, windowSize, windowOverlap, sequences, pilesPos[jobsLoaded].first, pilesPos[jobsLoaded].second, merSize, maxSupport, commonKMers);
		templates[jobsLoaded] = curPile[0];
    	results[jobsLoaded] = myPool.push(computeConsensusAssemblyPolishing, readId, curPile, pilesPos[jobsLoaded], minSupport, merSize, commonKMers, minAnchors, solidThresh, windowSize, maxMSA, path, nbThreads);
        jobsLoaded++;
	}

	// Load the remaining jobs as other jobs terminate
	int curJob = 0;
    std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>> curRes;
    while(jobsLoaded < jobsToProcess) {
    	// Get the job results
        curRes = results[curJob].get();
        consensuses[jobsCompleted] = curRes.first;
        merCounts[jobsCompleted] = curRes.second;
        jobsCompleted++;
        
        // Load the next job
        curPile = getAlignmentWindowsSequences(alignments, minSupport, windowSize, windowOverlap, sequences, pilesPos[jobsLoaded].first, pilesPos[jobsLoaded].second, merSize, maxSupport, commonKMers);
		templates[jobsLoaded] = curPile[0];
    	results[curJob] = myPool.push(computeConsensusAssemblyPolishing, readId, curPile, pilesPos[jobsLoaded], minSupport, merSize, commonKMers, minAnchors, solidThresh, windowSize, maxMSA, path, nbThreads);
        jobsLoaded++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

	// Wait for the remaining jobs to terminate
	while(jobsCompleted < jobsLoaded) {
        // Get the job results
        curRes = results[curJob].get();
        consensuses[jobsCompleted] = curRes.first;
        merCounts[jobsCompleted] = curRes.second;
        jobsCompleted++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

	// Align computed consensuses to the read
	std::string correctedRead = alignConsensus(readId, sequences[alignments[0].qName], consensuses, merCounts, pilesPos, templates, pilesPos[0].first, windowSize, windowOverlap, solidThresh, merSize);

	// Trim read if needed (ie when performing correction), and drop it if it contains too many uncorrected bases
	if (doTrimRead) {
		correctedRead = trimRead(correctedRead, 1);
		if (!dropRead(correctedRead)) {
			return std::make_pair(readId, correctedRead);
		} else {
			return std::make_pair(readId, "");
		}
	} else {
		return std::make_pair(readId, correctedRead);
	}
}

void runCorrection(std::string PAFIndex, std::string alignmentFile, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile, std::string proofFile, unsigned maxMSA, std::string path) {
	std::ifstream templates(PAFIndex);
	std::ifstream alignments(alignmentFile);
	std::vector<Overlap> curReadAlignments;
	std::string curRead, line;
	curRead = "";

	indexReads(readIndex, readsFile);
	if (proofFile != "") {
		indexReads(readIndex, proofFile)
;	}

	std::string curTpl;
	std::pair<std::string, std::string> curRes;

	while (!alignments.eof()) {
		curReadAlignments = getNextReadPile(alignments, maxSupport);
        while (curReadAlignments.size() == 0 and !alignments.eof()) {
        	curReadAlignments = getNextReadPile(alignments, maxSupport);
        }

        if (curReadAlignments.size() != 0) {
        	curRes = processContig(curReadAlignments, minSupport, maxSupport, windowSize, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, maxMSA, path, nbThreads);
		    if (curRes.second.length() != 0) {
		        std::cout << ">" << curRes.first << std::endl << curRes.second << std::endl;
		    }
        }
	}

}
