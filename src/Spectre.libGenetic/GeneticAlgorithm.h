#pragma once
#include "DataTypes.h"
#include "Generation.h"
#include "MutationOperator.h"
#include "CrossoverOperator.h"
#include "Classifier.h"
#include "Scorer.h"

namespace Spectre::libGenetic
{
class GeneticAlgorithm
{
public:
	GeneticAlgorithm(const Dataset* data, MutationOperator mutation, CrossoverOperator crossover, Selection selection, Scorer classifier, long generationSize);
	virtual ~GeneticAlgorithm() = default;

private:
	const Dataset* m_Data;
	MutationOperator m_Mutation;
	CrossoverOperator m_Crossover;
	Selection m_Selection;
	Scorer m_Scorer;
	Generation m_CurrentGeneration, m_NewGeneration;
	long m_GenerationSize;
};
}
