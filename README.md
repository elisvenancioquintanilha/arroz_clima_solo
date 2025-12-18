# Projeto: Integração de Dados Climáticos e de Solo para Cultivo de Arroz

## 1. Objetivo

Este projeto tem como objetivo integrar dados de clima e solo para análises de desempenho de cultivares de arroz no Brasil. A base final combina informações de:

- Condições climáticas por fases do ciclo do arroz
- Características do solo (AD_UM, Classe_AD, etc.)
- Dados experimentais de ensaios de arroz (genótipos, localização, repetições, etc.)

O pipeline inclui aquisição de dados climáticos da NASA POWER, processamento de dados de solo a partir de shapefiles do AD_Brasil, e junção final com os dados da embrapa.

---

## 2. Estrutura do Repositório

scripts/
  01_climate_by_stage.R      # Extração e estatísticas climáticas
  02_solo_processing.R       # Processamento dos dados de solo
  03_merge_arroz_solo_clima.R # Junção dos dados
data/
  dados_arroz.csv




---

## 3. Scripts

### 3.1 `01_climate_by_stage.R`
- Objetivo: extrair dados climáticos por fase do ciclo do arroz (vegetativa, reprodutiva e enchimento de grão) e calcular estatísticas agregadas.
- Fontes: NASA POWER API (`nasapower`).
- Principais etapas:
  1. Carregamento dos dados experimentais (`dados_arroz.csv`).
  2. Extração de coordenadas únicas e datas do ciclo fenológico.
  3. Funções para:
     - Baixar dados climáticos diários (`get_climate_data`)
     - Calcular estatísticas por fase (`calc_stats`)
     - Processar cada linha do dataset (`process_climate_stats`)
  4. Processamento paralelo com `furrr` para otimização.
  5. Resultado: arquivo `climate_data_dezembro.csv`.

### 3.2 `02_solo_processing.R`
- Objetivo: processar os dados de solo a partir de shapefiles do AD_Brasil, ajustando pontos que caem em "Área Edificada" e preparando os dados para junção com dados climáticos.
- Principais etapas:
  1. Carregamento e filtragem do shapefile do Brasil.
  2. Interseção com estados de interesse (Mato Grosso, Goiás, Pará, Rondônia, Piauí, Tocantins, Maranhão).
  3. Transformação em `sf` para compatibilidade com `ggplot2` e manipulação espacial.
  4. Ajuste de coordenadas problemáticas (perturbação de pontos para evitar "Área Edificada").
  5. Criação de colunas de resumo para classificação de solos (`Simb`, `AD_UM_cat`).
  6. Exportação: `dados_clima_solo_v03_dezembro.xlsx` e `dados_final_corrigido_completo.shp`.

### 3.3 `03_merge_arroz_solo_clima.R`
- Objetivo: integrar os dados experimentais de arroz com os dados climáticos e de solo processados.
- Principais etapas:
  1. Leitura de `dados_arroz.csv` e aplicação da função de pré-processamento (`data_processing_elis2.R`).
  2. Leitura de `dados_clima_solo_v03_dezembro.xlsx`.
  3. Ajuste dos formatos de data e nomes de colunas.
  4. Junção dos datasets via `left_join` usando chaves: `LATITUDE`, `LONGITUDE`, `DATE`, `DTF_data`, `PI_data`, `PM_data`.
  5. Reorganização das colunas de solo (`AD_UM`, `Classe_AD`, `Simb`, `AD_UM_cat`).
  6. Limpeza final:
     - Remoção de variáveis com muitos `NA`.
     - Remoção de variáveis sem sentido ou não climáticas.
     - Filtragem de observações com erros ou inconsistentes.
  7. Resultado final:
     - Arquivo Excel: `dados_completos_dezembro_bruto.xlsx`
     - Base final tratada: `dados_dezembro.xlsx`

---

## 4. Pacotes utilizados

- `tidyverse`, `dplyr`, `readxl`, `writexl`, `stringr`
- `nasapower` (extração de dados climáticos)
- `lubridate` (manipulação de dados)
- `furrr`, `future` (processamento paralelo)
- `sf`, `terra`, `geobr`, `ggplot2`, `purrr` (manipulação e visualização)

---

## 5. Observações importantes

1. Os dados climáticos são diários e foram agregados por fase do ciclo do arroz, com estatísticas como:
   - Média, mediana, IQR, cumulativa, proporção acima do percentil 90 e abaixo do 10.
2. Ajustes de coordenadas foram aplicados para evitar pontos em "Área Edificada".
3. As categorias de AD_UM foram classificadas via k-means em três níveis: Baixa, Média e Alta.
4. Todo o pipeline é reprodutível e modular, permitindo atualização com novos dados ou diferentes regiões.
5. O README serve como guia geral; comentários nos scripts explicam decisões específicas de processamento.

---

## 6. Estrutura dos dados finais

- `GEN`: Genótipo do arroz
- `REP`: Repetição do ensaio
- `TRIAL`: Nome do ensaio
- `LATITUDE` / `LONGITUDE`: Coordenadas do ensaio
- `DATE`, `PI_data`, `DTF_data`, `PM_data`: Datas do ciclo fenológico
- Variáveis climáticas agregadas por fase (`veg_*`, `repro_*`, `gf_*`)
- Variáveis de solo (`AD_UM`, `Classe_AD`, `AD_UM_cat`, `Simb`)
- Outras variáveis relevantes (ex: `GY`, `PHT`)

---

## 7. Contato
- E-mail: [elisvenancio22@gmail.com]


