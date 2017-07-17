#include "Generation.h"
#include "GeneticAlgorithm.h"

namespace Spectre::libGenetic
{
GeneticAlgorithm::GeneticAlgorithm(OffspringGenerator&& offspringGenerator, Scorer&& scorer, StopCondition&& stopCondition)
    : m_OffspringGenerator(offspringGenerator),
      m_Scorer(scorer),
      m_StopCondition(stopCondition)
{
}

Generation GeneticAlgorithm::evolve(Generation&& generation)
{
    while(!m_StopCondition())
    {
        const auto scores = m_Scorer.Score(generation);
        generation = m_OffspringGenerator.next(generation, scores);
    }
    return generation;
}
}
