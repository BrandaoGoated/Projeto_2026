# Projeto 2026 — Deteção de enzimas degradadoras de PET em metagenomas

Projeto da unidade curricular **Projeto em Bioinformática** (Mestrado em Bioinformática, Universidade do Minho).

Aplicação e validação da ferramenta [M-PARTY](https://github.com/ozefreitas/M-PARTY) para detetar enzimas degradadoras de **PET** em metagenomas ambientais, usando enzimas de referência da base de dados **PAZy**.

**Autor:** David Brandão
**Orientadores:** José Freitas, Diogo Cachetas, Andreia Salvador

## O que foi feito

A partir de enzimas PET-ativas do PAZy, obtiveram-se as sequências codificantes, que foram usadas como referência para procurar enzimas semelhantes em dois metagenomas (BioProject PRJNA849162) por mapeamento de reads com KMA. Os candidatos foram validados por BLAST e AlphaFold.

Foram detetados dois candidatos: um semelhante à PETase **PD3** (*Brucella anthropi*) e outro à esterase **PpEst** (*Ectopseudomonas oleovorans*), ambos confirmados por BLAST.

## Estrutura

```
Fase 1/    — Excel com todos os dados utilizados para o primeiro artigo
Fase2e3_ResultadosEFicheiros/   — todos os ficheiros relacionados ao M-Party
Fase2e3_ResultadosEFicheiros/Resources     — Principais resultados obtidos (contudo estes não contém os FASTQ visto que são 150GB de ficheiros)
```

## Ferramentas

M-PARTY · BLAST · AlphaFold
