# SIDLE, una aproximación a caracterizar la microbiota de muestras parafinadas (FFPE) 
Código usado para la realización del TFM. Se usó docker para crear un contenedor estable con las versiones compatibles de todos los programas y se usó la plataforma bioinoformática de análisis de microbiota QIIME2.

docker exec -it charming_mcclintock bash
conda activate qiime2-amplicon-2024.10



# Mama de Nejman
# 1. Preparación pre-sidle:
# 1.1: Creación de manifest.tsv y copiado a la máquina virtual (en PowerShell):
docker cp C:\Users\Sergio\Documents\dockermicrobiota\datanejman\Mama\Muestras charming_mcclintock:/data/nejman/Mama/
docker cp C:\Users\Sergio\Documents\dockermicrobiota\datanejman\archivos_mama_020625\Mama\metadata_mama.tsv charming_mcclintock:/data/nejman/Mama/
docker cp C:\Users\Sergio\Documents\dockermicrobiota\datanejman\Mama\manifest_mama.csv charming_mcclintock:/data/nejman/Mama/



# 1.2: Cambiar el nombre de los archivos forward y reverse con su id:
for dir in SRR*; do
  if [ -d "$dir" ]; then
    mv "$dir/forward.fastqsanger.gz" "$dir/${dir}_forward.fastq.gz"
    mv "$dir/reverse.fastqsanger.gz" "$dir/${dir}_reverse.fastq.gz"
  fi
done



# 1.3: Lanzar fastqc a las muestras:
fastqc SRR115*/*.gz -o fastqc -t 4



# 1.4: Lanzar multiqc a las muestras forward y reverse:
multiqc fastqc/SRR*_forward* -o multiqc/multiqc_f
multiqc fastqc/SRR*_reverse* -o multiqc/multiqc_r
multiqc fastqc/SRR* -o multiqc
  # Hay que quitar los adaptadores con cutadapt. 



# 1.5: Lanzar cutadapt para eliminar los adaptadores
for dir in SRR*/; do
  id=$(basename "$dir")
  r1="${dir}${id}_forward.fastq.gz"
  r2="${dir}${id}_reverse.fastq.gz"
  fixed_r1="${dir}${id}_fixed_forward.fastq"
  fixed_r2="${dir}${id}_fixed_reverse.fastq"
  out1="${dir}${id}_trimmed_forward.fastq.gz"
  out2="${dir}${id}_trimmed_reverse.fastq.gz"

  # Arreglar encabezados
  zcat "$r1" | sed 's/_forward\/1$/\/1/' > "$fixed_r1"
  zcat "$r2" | sed 's/_reverse\/2$/\/2/' > "$fixed_r2"

  # Ejecutar cutadapt
  cutadapt \
    -a AGATCGGAAGAG \
    -A AGATCGGAAGAG \
    -o "$out1" \
    -p "$out2" \
    "$fixed_r1" "$fixed_r2"

  # Eliminar archivos temporales
  rm "$fixed_r1" "$fixed_r2"
done

  #Eliminar muestras:
for dir in SRR*/; do
  find "$dir" -type f -name '*trimmed*' -exec rm -v {} \;
done



# 1.6: Lanzar fastqc a las muestras:
fastqc SRR*/*trimmed* -o fastqc_trimmed -t 4



# 1.7: Lanzar multiqc a las muestras forward y reverse:
multiqc fastqc_trimmed/SRR*_trimmed_forward* -o multiqc_trimmed/multiqcForward_trimmed --force
multiqc fastqc_trimmed/SRR*_trimmed_reverse* -o multiqc_trimmed/multiqcReverse_trimmed --force
multiqc fastqc_trimmed/SRR* -o multiqc_trimmed --force
  # Cortar en 124



# 1.8: Transformamos las muestras demultiplexeadas al archivo .qza
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33 \
  --input-path ./manifest_mama.csv \
  --output-path ./demux/demux-paired-end.qza
  ## En este caso el manifest es .csv y tiene las siguientes columnas: sample-id,	absolute-filepath,	direction

  ## Comprobación
qiime demux summarize \
  --i-data ./demux/demux-paired-end.qza \
  --o-visualization ./demux/demux-paired-end.qzv



# 1.9: Demultiplexear en regiones
# 1.9.1: Región FR1 (V2)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ./demux/demux-paired-end.qza \
 --p-front-f TGGCGAACGGGTGAGTAA \
 --p-front-r CCGTGTCTCAGTCCCARTG \
 --p-error-rate 0.1 \
 --p-indels \
 --p-discard-untrimmed \
 --o-trimmed-sequences demux/FR1/FR1_seqs.qza

    ## Comprobación
qiime demux summarize \
  --i-data demux/FR1/FR1_seqs.qza \
  --o-visualization demux/FR1/FR1_seqs.qzv


# 1.9.2: Región FR2 (V3)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ./demux/demux-paired-end.qza \
 --p-front-f ACTCCTACGGGAGGCAGC \
 --p-front-r GTATTACCGCGGCTGCTG \
 --p-error-rate 0.1 \
 --p-indels \
 --p-discard-untrimmed \
 --o-trimmed-sequences demux/FR2/FR2_seqs.qza
    
    ## Comprobación
qiime demux summarize \
  --i-data demux/FR2/FR2_seqs.qza \
  --o-visualization demux/FR2/FR2_seqs.qzv


# 1.9.3: Región FR3 (V5)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ./demux/demux-paired-end.qza \
 --p-front-f GTGTAGCGGTGRAATGCG \
 --p-front-r CCCGTCAATTCMTTTGAGTT \
 --p-error-rate 0.1 \
 --p-indels \
 --p-discard-untrimmed \
 --o-trimmed-sequences demux/FR3/FR3_seqs.qza
    
    ## Comprobación
qiime demux summarize \
  --i-data demux/FR3/FR3_seqs.qza \
  --o-visualization demux/FR3/FR3_seqs.qzv

  
# 1.9.4: Región FR4 (V6)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ./demux/demux-paired-end.qza \
 --p-front-f GGAGCATGTGGWTTAATTCGA \
 --p-front-r CGTTGCGGGACTTAACCC \
 --p-error-rate 0.1 \
 --p-indels \
 --p-discard-untrimmed \
 --o-trimmed-sequences demux/FR4/FR4_seqs.qza
    
    ## Comprobación
qiime demux summarize \
  --i-data demux/FR4/FR4_seqs.qza \
  --o-visualization demux/FR4/FR4_seqs.qzv

  
# 1.9.5: Región FR5 (V8)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences ./demux/demux-paired-end.qza \
 --p-front-f GGAGGAAGGTGGGGATGAC \
 --p-front-r AAGGCCCGGGAACGTATT \
 --p-error-rate 0.1 \
 --p-indels \
 --p-discard-untrimmed \
 --o-trimmed-sequences demux/FR5/FR5_seqs.qza
    
    ## Comprobación
qiime demux summarize \
  --i-data demux/FR5/FR5_seqs.qza \
  --o-visualization demux/FR5/FR5_seqs.qzv


# 1.10: Denoising con DADA2
# 1.10.1: FR1 (V2)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux/FR1/FR1_seqs.qza \
  --p-trunc-len-f 106 \
  --p-trunc-len-r 106 \
  --o-table demux/FR1/FR1_table.qza \
  --o-representative-sequences demux/FR1/FR1_rep_seqs.qza \
  --o-denoising-stats demux/FR1/FR1_stats.qza

  ## Visualizar tabla de frecuencias
qiime feature-table summarize \
  --i-table demux/FR1/FR1_table.qza \
  --o-visualization demux/FR1/FR1_table.qzv \
  --m-sample-metadata-file metadata_mama.tsv

  ## Visualizar secuencias representativas
qiime feature-table tabulate-seqs \
  --i-data demux/FR1/FR1_rep_seqs.qza \
  --o-visualization demux/FR1/FR1_rep_seqs.qzv


# 1.10.2: FR2 (V3)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux/FR2/FR2_seqs.qza \
  --p-trunc-len-f 106 \
  --p-trunc-len-r 106 \
  --o-table demux/FR2/FR2_table.qza \
  --o-representative-sequences demux/FR2/FR2_rep_seqs.qza \
  --o-denoising-stats demux/FR2/FR2_stats.qza

  ## Visualizar tabla de frecuencias
qiime feature-table summarize \
  --i-table demux/FR2/FR2_table.qza \
  --o-visualization demux/FR2/FR2_table.qzv \
  --m-sample-metadata-file metadata_mama.tsv

  ## Visualizar secuencias representativas
qiime feature-table tabulate-seqs \
  --i-data demux/FR2/FR2_rep_seqs.qza \
  --o-visualization demux/FR2/FR2_rep_seqs.qzv


# 1.10.3: FR3 (V5)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux/FR3/FR3_seqs.qza \
  --p-trim-left-f 6 \
  --p-trim-left-r 6 \
  --p-trunc-len-f 107 \
  --p-trunc-len-r 105 \
  --o-table demux/FR3/FR3_table.qza \
  --o-representative-sequences demux/FR3/FR3_rep_seqs.qza \
  --o-denoising-stats demux/FR3/FR3_stats.qza

  ## Visualizar tabla de frecuencias
qiime feature-table summarize \
  --i-table demux/FR3/FR3_table.qza \
  --o-visualization demux/FR3/FR3_table.qzv \
  --m-sample-metadata-file metadata_mama.tsv

  ## Visualizar secuencias representativas
qiime feature-table tabulate-seqs \
  --i-data demux/FR3/FR3_rep_seqs.qza \
  --o-visualization demux/FR3/FR3_rep_seqs.qzv


# 1.10.4: FR4 (V6)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux/FR4/FR4_seqs.qza \
  --p-trunc-len-f 104 \
  --p-trunc-len-r 107 \
  --o-table demux/FR4/FR4_table.qza \
  --o-representative-sequences demux/FR4/FR4_rep_seqs.qza \
  --o-denoising-stats demux/FR4/FR4_stats.qza

  ## Visualizar tabla de frecuencias
qiime feature-table summarize \
  --i-table demux/FR4/FR4_table.qza \
  --o-visualization demux/FR4/FR4_table.qzv \
  --m-sample-metadata-file metadata_mama.tsv

  ## Visualizar secuencias representativas
qiime feature-table tabulate-seqs \
  --i-data demux/FR4/FR4_rep_seqs.qza \
  --o-visualization demux/FR4/FR4_rep_seqs.qzv


# 1.10.2: FR5 (V8)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux/FR5/FR5_seqs.qza \
  --p-trunc-len-f 106 \
  --p-trunc-len-r 107 \
  --o-table demux/FR5/FR5_table.qza \
  --o-representative-sequences demux/FR5/FR5_rep_seqs.qza \
  --o-denoising-stats demux/FR5/FR5_stats.qza

  ## Visualizar tabla de frecuencias
qiime feature-table summarize \
  --i-table demux/FR5/FR5_table.qza \
  --o-visualization demux/FR5/FR5_table.qzv \
  --m-sample-metadata-file metadata_mama.tsv

  ## Visualizar secuencias representativas
qiime feature-table tabulate-seqs \
  --i-data demux/FR5/FR5_rep_seqs.qza \
  --o-visualization demux/FR5/FR5_rep_seqs.qzv



# 1.11: Recortar las secuencias representativas (sidle usa secuencias de misma longitud)
# 1.11.1: FR1 (V2)
qiime sidle trim-dada2-posthoc \
 --i-table demux/FR1/FR1_table.qza \
 --i-representative-sequences demux/FR1/FR1_rep_seqs.qza \
 --p-trim-length 193 \
 --o-trimmed-table demux/FR1/FR1_table_193.qza \
 --o-trimmed-representative-sequences demux/FR1/FR1_rep-seqs_193.qza
 # 193/195/196

  ## Visualizar secuencias representativas recortadas
qiime feature-table tabulate-seqs \
  --i-data demux/FR1/FR1_rep-seqs_193.qza \
  --o-visualization demux/FR1/FR1_rep-seqs_193.qzv


# 1.11.2: FR2 (V3)
qiime sidle trim-dada2-posthoc \
 --i-table demux/FR2/FR2_table.qza \
 --i-representative-sequences demux/FR2/FR2_rep_seqs.qza \
 --p-trim-length 142 \
 --o-trimmed-table demux/FR2/FR2_table_142.qza \
 --o-trimmed-representative-sequences demux/FR2/FR2_rep-seqs_142.qza
 # 142/158/163

  ## Visualizar secuencias representativas recortadas
qiime feature-table tabulate-seqs \
  --i-data demux/FR2/FR2_rep-seqs_142.qza \
  --o-visualization demux/FR2/FR2_rep-seqs_142.qzv


# 1.11.3: FR3 (V5)
qiime sidle trim-dada2-posthoc \
 --i-table demux/FR3/FR3_table.qza \
 --i-representative-sequences demux/FR3/FR3_rep_seqs.qza \
 --p-trim-length 188 \
 --o-trimmed-table demux/FR3/FR3_table_188.qza \
 --o-trimmed-representative-sequences demux/FR3/FR3_rep-seqs_188.qza
 # 188. Solo hay 10 secuencias

  ## Visualizar secuencias representativas recortadas
qiime feature-table tabulate-seqs \
  --i-data demux/FR3/FR3_rep-seqs_188.qza \
  --o-visualization demux/FR3/FR3_rep-seqs_188.qzv


# 1.11.4: FR4 (V6)
qiime sidle trim-dada2-posthoc \
 --i-table demux/FR4/FR4_table.qza \
 --i-representative-sequences demux/FR4/FR4_rep_seqs.qza \
 --p-trim-length 121 \
 --o-trimmed-table demux/FR4/FR4_table_121.qza \
 --o-trimmed-representative-sequences demux/FR4/FR4_rep-seqs_121.qza
 # 121/122/124

  ## Visualizar secuencias representativas recortadas
qiime feature-table tabulate-seqs \
  --i-data demux/FR4/FR4_rep-seqs_121.qza \
  --o-visualization demux/FR4/FR4_rep-seqs_121.qzv


# 1.11.5: FR5 (V8)
qiime sidle trim-dada2-posthoc \
 --i-table demux/FR5/FR5_table.qza \
 --i-representative-sequences demux/FR5/FR5_rep_seqs.qza \
 --p-trim-length 181 \
 --o-trimmed-table demux/FR5/FR5_table_181.qza \
 --o-trimmed-representative-sequences demux/FR5/FR5_rep-seqs_181.qza
 # 181/181/182

  ## Visualizar secuencias representativas recortadas
qiime feature-table tabulate-seqs \
  --i-data demux/FR5/FR5_rep-seqs_181.qza \
  --o-visualization demux/FR5/FR5_rep-seqs_181.qzv

-------------------------------------------------------------------------------------------------------------------------------------------

# 2. Preparación Sidle:
# 2.1: Importar la base de datos (GreenGenes 13.8) en Qiime2
# 2.1.1: Importar las secuencias
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path Greengenes/99_otus.fasta \
  --output-path Greengenes/99_otus.qza \
  --input-format DNAFASTAFormat


# 2.1.2: Importar la taxonomía 
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path Greengenes/99_otu_taxonomy.txt \
  --output-path Greengenes/99_otu_taxonomy.qza



# 2.2: Eliminar las secuencias con 3 o más nucleotidos degenerados de la base de datos
 qiime rescript cull-seqs \
 --i-sequences Greengenes/99_otus.qza \
 --p-num-degenerates 3 \
 --o-clean-sequences dataBase/filtered_db.qza

  ## Comprobación del filtrado
 qiime feature-table tabulate-seqs \
 --i-data dataBase/filtered_db.qza \
 --o-visualization dataBase/filtered_db.qzv



# 2.3: Preparación de la base de datos para los primers. Utilizaron el conjunto de primers 5R sin especificar regiones.
# 2.3.1: FR1 (V2)
  ## Extracción de la región FR1
qiime feature-classifier extract-reads \
 --i-sequences dataBase/filtered_db.qza \
 --p-f-primer TGGCGAACGGGTGAGTAA  \
 --p-r-primer CCGTGTCTCAGTCCCARTG  \
 --o-reads dataBase/FR1/FR1_filtered_db.qza

  ## Preparación de FR1 para Sidle
qiime sidle prepare-extracted-region \
  --i-sequences dataBase/FR1/FR1_filtered_db.qza \
  --p-region "FR1" \
  --p-fwd-primer TGGCGAACGGGTGAGTAA \
  --p-rev-primer CCGTGTCTCAGTCCCARTG  \
  --p-trim-length 193 \
  --o-collapsed-kmers dataBase/FR1/FR1_filtered_db_193_kmers.qza \
  --o-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza 

  ## Visualización del Kmer-map:
qiime metadata tabulate \
 --m-input-file dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
 --o-visualization dataBase/FR1/FR1_filtered_db_193_kmer-map.qzv


# 2.3.2: FR2 (V3)
  ## Extracción de la región FR2
qiime feature-classifier extract-reads \
 --i-sequences dataBase/filtered_db.qza \
 --p-f-primer ACTCCTACGGGAGGCAGC  \
 --p-r-primer GTATTACCGCGGCTGCTG \
 --o-reads dataBase/FR2/FR2_filtered_db.qza

  ## Preparación de FR2 para Sidle
qiime sidle prepare-extracted-region \
  --i-sequences dataBase/FR2/FR2_filtered_db.qza \
  --p-region "FR2" \
  --p-fwd-primer ACTCCTACGGGAGGCAGC \
  --p-rev-primer GTATTACCGCGGCTGCTG \
  --p-trim-length 142 \
  --o-collapsed-kmers dataBase/FR2/FR2_filtered_db_142_kmers.qza \
  --o-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza 

  ## Visualización del Kmer-map:
qiime metadata tabulate \
 --m-input-file dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
 --o-visualization dataBase/FR2/FR2_filtered_db_142_kmer-map.qzv


# 2.3.3: FR3 (V5)
  ## Extracción de la región FR3
qiime feature-classifier extract-reads \
 --i-sequences dataBase/filtered_db.qza \
 --p-f-primer GTGTAGCGGTGRAATGCG  \
 --p-r-primer CCCGTCAATTCMTTTGAGTT \
 --o-reads dataBase/FR3/FR3_filtered_db.qza

  ## Preparación de FR3 para Sidle
qiime sidle prepare-extracted-region \
  --i-sequences dataBase/FR3/FR3_filtered_db.qza \
  --p-region "FR3" \
  --p-fwd-primer GTGTAGCGGTGRAATGCG \
  --p-rev-primer CCCGTCAATTCMTTTGAGTT \
  --p-trim-length 188 \
  --o-collapsed-kmers dataBase/FR3/FR3_filtered_db_188_kmers.qza \
  --o-kmer-map dataBase/FR3/FR3_filtered_db_188_kmer-map.qza 

  ## Visualización del Kmer-map:
qiime metadata tabulate \
 --m-input-file dataBase/FR3/FR3_filtered_db_188_kmer-map.qza  \
 --o-visualization dataBase/FR3/FR3_filtered_db_188_kmer-map.qzV


 # 2.3.4: FR4 (V6)
  ## Extracción de la región FR4
qiime feature-classifier extract-reads \
 --i-sequences dataBase/filtered_db.qza \
 --p-f-primer GGAGCATGTGGWTTAATTCGA  \
 --p-r-primer CGTTGCGGGACTTAACCC \
 --o-reads dataBase/FR4/FR4_filtered_db.qza

  ## Preparación de FR4 para Sidle
qiime sidle prepare-extracted-region \
  --i-sequences dataBase/FR4/FR4_filtered_db.qza \
  --p-region "FR4" \
  --p-fwd-primer GGAGCATGTGGWTTAATTCGA \
  --p-rev-primer CGTTGCGGGACTTAACCC \
  --p-trim-length 121 \
  --o-collapsed-kmers dataBase/FR4/FR4_filtered_db_121_kmers.qza \
  --o-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza 

  ## Visualización del Kmer-map:
qiime metadata tabulate \
 --m-input-file dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
 --o-visualization dataBase/FR4/FR4_filtered_db_121_kmer-map.qzV


 # 2.3.5: FR5 (V8)
  ## Extracción de la región FR5
qiime feature-classifier extract-reads \
 --i-sequences dataBase/filtered_db.qza \
 --p-f-primer GGAGGAAGGTGGGGATGAC  \
 --p-r-primer AAGGCCCGGGAACGTATT \
 --o-reads dataBase/FR5/FR5_filtered_db.qza

  ## Preparación de FR5 para Sidle
qiime sidle prepare-extracted-region \
  --i-sequences dataBase/FR5/FR5_filtered_db.qza \
  --p-region "FR5" \
  --p-fwd-primer GGAGGAAGGTGGGGATGAC \
  --p-rev-primer AAGGCCCGGGAACGTATT \
  --p-trim-length 181 \
  --o-collapsed-kmers dataBase/FR5/FR5_filtered_db_181_kmers.qza \
  --o-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza 

  ## Visualización del Kmer-map:
qiime metadata tabulate \
 --m-input-file dataBase/FR5/FR5_filtered_db_181_kmer-map.qza  \
 --o-visualization dataBase/FR5/FR5_filtered_db_181_kmer-map.qzv


# -------------------------------------------------------------------------------------------------------------------------------------------

# 3. Reconstrucción Sidle:
# 3.1: Alineamiento Regional
  ## FR1
 qiime sidle align-regional-kmers \
--i-kmers dataBase/FR1/FR1_filtered_db_193_kmers.qza \
--i-rep-seq demux/FR1/FR1_rep-seqs_193.qza \
--p-region FR1 \
--p-n-workers 2 \
--o-regional-alignment sidle/alignment/FR1_193_align-map.qza


  ## FR2
 qiime sidle align-regional-kmers \
--i-kmers dataBase/FR2/FR2_filtered_db_142_kmers.qza \
--i-rep-seq demux/FR2/FR2_rep-seqs_142.qza \
--p-region FR2 \
--p-n-workers 2 \
--o-regional-alignment sidle/alignment/FR2_142_align-map.qza


  ## FR3
 qiime sidle align-regional-kmers \
--i-kmers dataBase/FR3/FR3_filtered_db_188_kmers.qza \
--i-rep-seq demux/FR3/FR3_rep-seqs_188.qza \
--p-region "FR3" \
--p-n-workers 2 \
--o-regional-alignment sidle/alignment/FR3_188_align-map.qza


  ## FR4
 qiime sidle align-regional-kmers \
--i-kmers dataBase/FR4/FR4_filtered_db_121_kmers.qza \
--i-rep-seq demux/FR4/FR4_rep-seqs_121.qza \
--p-region FR4 \
--p-n-workers 2 \
--o-regional-alignment sidle/alignment/FR4_121_align-map.qza


  ## FR5
 qiime sidle align-regional-kmers \
--i-kmers dataBase/FR5/FR5_filtered_db_181_kmers.qza \
--i-rep-seq demux/FR5/FR5_rep-seqs_181.qza \
--p-region FR5 \
--p-n-workers 2 \
--o-regional-alignment sidle/alignment/FR5_181_align-map.qza

qiime tools export --input-path sidle/alignment/FR5_181_align-map.qza --output-path export_FR5


# 3.2: Reconstrucción
# 3.2.1: Reconstrucción de la base de datos
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR3 \
  --i-kmer-map dataBase/FR3/FR3_filtered_db_188_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR3_188_align-map.qza \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/mama_recon_summary.qzv


# 3.2.2: Reconstrucción de la tabla
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR3 \
    --i-regional-alignment sidle/alignment/FR3_188_align-map.qza \
    --i-regional-table demux/FR3/FR3_table_188.qza \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/tabla_recon.qza \
 --o-visualization sidle/reconstruction/tabla_recon.qzv

# 3.2.3: Reconstrucción taxonómica
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/taxa-bar-plots.qzv


# -----------------------------------------------------------------------------------------------------------------------------------------------

# R4: Reconstrucción FR1/FR2/FR4/FR5 (RECONSTRUYE)
# R4.1: Reconstrucción de la base de datos FR1/FR2/FR4/FR5
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_1.2.4.5_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_1.2.4.5_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2.4.5_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_1.2.4.5_mama_recon_summary.qzv


# R4.2: Reconstrucción de la tabla FR1/FR2/FR4/FR5
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/FR_1.2.4.5_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_1.2.4.5_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_1.2.4.5_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_1.2.4.5_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_1.2.4.5_tabla_recon.qzv


# R4.3: Reconstrucción taxonómica FR1/FR2/FR4/FR5
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_1.2.4.5_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_1.2.4.5_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2.4.5_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_1.2.4.5_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_1.2.4.5_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_1.2.4.5_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_1.2.4.5_taxa-bar-plots.qzv


# R4.D: Análisis de diversidad
  ## Valorar el sampling-depth
qiime diversity alpha-rarefaction \
--i-table sidle/reconstruction/FR_1.2.4.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4.5/FR_1.2.4.5_rarefaction_curves.qzv \
--p-min-depth 1000 \
--p-max-depth 125694

  ## Core-metrics
qiime diversity core-metrics \
--i-table sidle/reconstruction/FR_1.2.4.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--p-sampling-depth 28000 \
--output-dir sidle/reconstruction/diversidad/FR_1.2.4.5

  ## Diversidad Alfa
qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.2.4.5/observed_features_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4.5/observed_features_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.2.4.5/shannon_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4.5/shannon_vector.qzv

  ## Diversidad Beta
qiime diversity beta-group-significance \
--i-distance-matrix sidle/reconstruction/diversidad/FR_1.2.4.5/bray_curtis_distance_matrix.qza \
--m-metadata-file metadata_mama.tsv \
--m-metadata-column Subtype \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4.5/bray_curtis_distance_matrix_Subtype.qzv \
--p-pairwise 

# -----------------------------------------------------------------------------------------------------------------------------------------------

# R2: Reconstrucción FR1/FR2 (NO RECONSTRUYE)
# R2.1: Reconstrucción de la base de datos FR1/FR2
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_1.2_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_1.2_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_1.2_mama_recon_summary.qzv


# R2.2: Reconstrucción de la tabla FR1/FR2
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --i-database-map sidle/reconstruction/FR_1.2_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_1.2_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_1.2_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_1.2_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_1.2_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR1/FR2
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_1.2_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_1.2_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_1.2_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_1.2_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_1.2_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_1.2_taxa-bar-plots.qzv


# -----------------------------------------------------------------------------------------------------------------------------------------------

# R2: Reconstrucción FR1/FR4 (RECONSTRUYE)
# R2.1: Reconstrucción de la base de datos FR1/FR4
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_1.4_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_1.4_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.4_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_1.4_mama_recon_summary.qzv


# R2.2: Reconstrucción de la tabla FR1/FR4
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --i-database-map sidle/reconstruction/FR_1.4_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_1.4_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_1.4_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_1.4_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_1.4_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR1/FR4
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_1.4_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_1.4_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.4_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_1.4_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_1.4_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_1.4_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_1.4_taxa-bar-plots.qzv


# R2.D: Análisis de diversidad
  ## Valorar el sampling-depth
qiime diversity alpha-rarefaction \
--i-table sidle/reconstruction/FR_1.4_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.4_rarefaction_curves.qzv \
--p-min-depth 900 \
--p-max-depth 129836

  ## Core-metrics
qiime diversity core-metrics \
--i-table sidle/reconstruction/FR_1.4_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--p-sampling-depth 11196 \
--output-dir sidle/reconstruction/diversidad/FR_1.4

  ## Diversidad Alfa
qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.4/observed_features_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.4/observed_features_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.4/shannon_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.4/shannon_vector.qzv

  ## Diversidad Beta
qiime diversity beta-group-significance \
--i-distance-matrix sidle/reconstruction/diversidad/FR_1.4/bray_curtis_distance_matrix.qza \
--m-metadata-file metadata_mama.tsv \
--m-metadata-column Subtype \
--o-visualization sidle/reconstruction/diversidad/FR_1.4/bray_curtis_distance_matrix_Subtype.qzv \
--p-pairwise 


  # -----------------------------------------------------------------------------------------------------------------------------------------------

# R2: Reconstrucción FR1/FR5 (NO RECONSTRUYE)
# R2.1: Reconstrucción de la base de datos FR1/FR5
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_1.5_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_1.5_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.5_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_1.5_mama_recon_summary.qzv


# R2.2: Reconstrucción de la tabla FR1/FR5
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/FR_1.5_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_1.5_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_1.5_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_1.5_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_1.5_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR1/FR5
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_1.5_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_1.5_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.5_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_1.5_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_1.5_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_1.5_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_1.5_taxa-bar-plots.qzv


 # -----------------------------------------------------------------------------------------------------------------------------------------------

  # R2: Reconstrucción FR2/FR4 (NO RECONSTRUYE)
# R2.1: Reconstrucción de la base de datos FR2/FR4
 qiime sidle reconstruct-database \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_2.4_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_2.4_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_2.4_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_2.4_mama_recon_summary.qzv


# R2.2: Reconstrucción de la tabla FR2/FR4
 qiime sidle reconstruct-counts \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --i-database-map sidle/reconstruction/FR_2.4_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_2.4_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_2.4_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_2.4_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_2.4_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR2/FR4
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_2.4_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_2.4_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_2.4_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_2.4_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_2.4_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_2.4_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_2.4_taxa-bar-plots.qzv


# -----------------------------------------------------------------------------------------------------------------------------------------------

  # R2: Reconstrucción FR2/FR5 (NO RECONSTRUYE)
# R2.1: Reconstrucción de la base de datos FR2/FR5
 qiime sidle reconstruct-database \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_2.5_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_2.5_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_2.5_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_2.5_mama_recon_summary.qzv


# R2.2: Reconstrucción de la tabla FR2/FR5
 qiime sidle reconstruct-counts \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/FR_2.5_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_2.5_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_2.5_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_2.5_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_2.5_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR2/FR5
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_2.5_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_2.5_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_2.5_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_2.5_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_2.5_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_2.5_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_2.5_taxa-bar-plots.qzv


# -----------------------------------------------------------------------------------------------------------------------------------------------

  # R2: Reconstrucción FR4/FR5 (RECONSTRUYE)
# R2.1: Reconstrucción de la base de datos FR4/FR5
 qiime sidle reconstruct-database \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_4.5_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_4.5_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_4.5_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_4.5_mama_recon_summary.qzv


# R2.2: Reconstrucción de la tabla FR4/FR5
 qiime sidle reconstruct-counts \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/FR_4.5_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_4.5_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_4.5_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_4.5_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_4.5_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR4/FR5
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_4.5_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_4.5_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_4.5_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_4.5_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_4.5_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_4.5_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_4.5_taxa-bar-plots.qzv


# R2.D: Análisis de diversidad
  ## Valorar el sampling-depth
qiime diversity alpha-rarefaction \
--i-table sidle/reconstruction/FR_4.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_4.5_rarefaction_curves.qzv \
--p-min-depth 900 \
--p-max-depth 191542

  ## Core-metrics
qiime diversity core-metrics \
--i-table sidle/reconstruction/FR_4.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--p-sampling-depth 21340 \
--output-dir sidle/reconstruction/diversidad/FR_4.5

  ## Diversidad Alfa
qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_4.5/observed_features_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_4.5/observed_features_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_4.5/shannon_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_4.5/shannon_vector.qzv

  ## Diversidad Beta
qiime diversity beta-group-significance \
--i-distance-matrix sidle/reconstruction/diversidad/FR_4.5/bray_curtis_distance_matrix.qza \
--m-metadata-file metadata_mama.tsv \
--m-metadata-column Subtype \
--o-visualization sidle/reconstruction/diversidad/FR_4.5/bray_curtis_distance_matrix_Subtype.qzv \
--p-pairwise 

# -----------------------------------------------------------------------------------------------------------------------------------------------

  # R3: Reconstrucción FR1/FR2/FR4 (RECONSTRUYE)
# R3.1: Reconstrucción de la base de datos FR1/FR2/FR4
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_1.2.4_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_1.2.4_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2.4_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_1.2.4_mama_recon_summary.qzv


# R3.2: Reconstrucción de la tabla FR1/FR2/FR4
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --i-database-map sidle/reconstruction/FR_1.2.4_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_1.2.4_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_1.2.4_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_1.2.4_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_1.2.4_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR1/FR2/FR4
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_1.2.4_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_1.2.4_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2.4_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_1.2.4_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_1.2.4_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_1.2.4_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_1.2.4_taxa-bar-plots.qzv


# R3.D: Análisis de diversidad
  ## Valorar el sampling-depth
qiime diversity alpha-rarefaction \
--i-table sidle/reconstruction/FR_1.2.4_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4_rarefaction_curves.qzv \
--p-min-depth 900 \
--p-max-depth 132870

  ## Core-metrics
qiime diversity core-metrics \
--i-table sidle/reconstruction/FR_1.2.4_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--p-sampling-depth 15016 \
--output-dir sidle/reconstruction/diversidad/FR_1.2.4

  ## Diversidad Alfa
qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.2.4/observed_features_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4/observed_features_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.2.4/shannon_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4/shannon_vector.qzv

  ## Diversidad Beta
qiime diversity beta-group-significance \
--i-distance-matrix sidle/reconstruction/diversidad/FR_1.2.4/bray_curtis_distance_matrix.qza \
--m-metadata-file metadata_mama.tsv \
--m-metadata-column Subtype \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.4/bray_curtis_distance_matrix_Subtype.qzv \
--p-pairwise 

# -----------------------------------------------------------------------------------------------------------------------------------------------

  # R3: Reconstrucción FR1/FR2/FR5 (RECONSTRUYE)
# R3.1: Reconstrucción de la base de datos FR1/FR2/FR5
 qiime sidle reconstruct-database \
 --p-region FR1 \
  --i-kmer-map dataBase/FR1/FR1_filtered_db_193_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 2 \
 --o-database-map sidle/reconstruction/FR_1.2.5_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_1.2.5_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2.5_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_1.2.5_mama_recon_summary.qzv


# R3.2: Reconstrucción de la tabla FR1/FR2/FR5
 qiime sidle reconstruct-counts \
  --p-region FR1 \
    --i-regional-alignment sidle/alignment/FR1_193_align-map.qza \
    --i-regional-table demux/FR1/FR1_table_193.qza \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/FR_1.2.5_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_1.2.5_mama_recon_summary.qza \
  --p-n-workers 2 \
  --o-reconstructed-table sidle/reconstruction/FR_1.2.5_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_1.2.5_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_1.2.5_tabla_recon.qzv


# R2.3: Reconstrucción taxonómica FR1/FR2/FR5
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_1.2.5_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_1.2.5_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_1.2.5_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_1.2.5_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_1.2.5_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_1.2.5_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_1.2.5_taxa-bar-plots.qzv


# R3.D: Análisis de diversidad
  ## Valorar el sampling-depth
qiime diversity alpha-rarefaction \
--i-table sidle/reconstruction/FR_1.2.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.5_rarefaction_curves.qzv \
--p-min-depth 300 \
--p-max-depth 136331

  ## Core-metrics
qiime diversity core-metrics \
--i-table sidle/reconstruction/FR_1.2.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--p-sampling-depth 18024 \
--output-dir sidle/reconstruction/diversidad/FR_1.2.5

  ## Diversidad Alfa
qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.2.5/observed_features_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.5/observed_features_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_1.2.5/shannon_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.5/shannon_vector.qzv

  ## Diversidad Beta
qiime diversity beta-group-significance \
--i-distance-matrix sidle/reconstruction/diversidad/FR_1.2.5/bray_curtis_distance_matrix.qza \
--m-metadata-file metadata_mama.tsv \
--m-metadata-column Subtype \
--o-visualization sidle/reconstruction/diversidad/FR_1.2.5/bray_curtis_distance_matrix_Subtype.qzv \
--p-pairwise 

# -----------------------------------------------------------------------------------------------------------------------------------------------

  # R3: Reconstrucción FR2/FR4/FR5 (RECONSTRUYE)
# R3.1: Reconstrucción de la base de datos FR2/FR4/FR5
 qiime sidle reconstruct-database \
 --p-region FR2 \
  --i-kmer-map dataBase/FR2/FR2_filtered_db_142_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
 --p-region FR4 \
  --i-kmer-map dataBase/FR4/FR4_filtered_db_121_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
 --p-region FR5 \
  --i-kmer-map dataBase/FR5/FR5_filtered_db_181_kmer-map.qza \
  --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
 --p-n-workers 1 \
 --o-database-map sidle/reconstruction/FR_2.4.5_mama_recon_map.qza \
 --o-database-summary sidle/reconstruction/FR_2.4.5_mama_recon_summary.qza

  ## Visualizar el summary para evaluar la calidad de la reconstrucción
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_2.4.5_mama_recon_summary.qza \
 --o-visualization sidle/reconstruction/FR_2.4.5_mama_recon_summary.qzv


# R3.2: Reconstrucción de la tabla FR2/FR4/FR5
 qiime sidle reconstruct-counts \
  --p-region FR2 \
    --i-regional-alignment sidle/alignment/FR2_142_align-map.qza \
    --i-regional-table demux/FR2/FR2_table_142.qza \
  --p-region FR4 \
    --i-regional-alignment sidle/alignment/FR4_121_align-map.qza \
    --i-regional-table demux/FR4/FR4_table_121.qza \
  --p-region FR5 \
    --i-regional-alignment sidle/alignment/FR5_181_align-map.qza \
    --i-regional-table demux/FR5/FR5_table_181.qza \
  --i-database-map sidle/reconstruction/FR_2.4.5_mama_recon_map.qza \
  --i-database-summary sidle/reconstruction/FR_2.4.5_mama_recon_summary.qza \
  --p-n-workers 1 \
  --o-reconstructed-table sidle/reconstruction/FR_2.4.5_tabla_recon.qza

  ## Visualizar el summary
qiime feature-table summarize \
 --i-table sidle/reconstruction/FR_2.4.5_tabla_recon.qza \
 --o-visualization sidle/reconstruction/FR_2.4.5_tabla_recon.qzv

# R2.3: Reconstrucción taxonómica FR2/FR4/FR5
qiime sidle reconstruct-taxonomy \
 --i-reconstruction-map sidle/reconstruction/FR_2.4.5_mama_recon_map.qza \
 --i-taxonomy Greengenes/99_otu_taxonomy.qza \
 --p-database 'greengenes' \
 --p-define-missing 'inherit' \
 --o-reconstructed-taxonomy sidle/reconstruction/FR_2.4.5_mama_recon_taxonomy.qza

  ## Visualizar la taxonomía
 qiime metadata tabulate \
 --m-input-file sidle/reconstruction/FR_2.4.5_mama_recon_taxonomy.qza \
 --o-visualization sidle/reconstruction/FR_2.4.5_mama_recon_taxonomy.qzv

 qiime taxa barplot \
  --i-table sidle/reconstruction/FR_2.4.5_tabla_recon.qza \
  --i-taxonomy sidle/reconstruction/FR_2.4.5_mama_recon_taxonomy.qza \
  --m-metadata-file metadata_mama.tsv \
  --o-visualization sidle/reconstruction/FR_2.4.5_taxa-bar-plots.qzv


# R3.D: Análisis de diversidad
  ## Valorar el sampling-depth
qiime diversity alpha-rarefaction \
--i-table sidle/reconstruction/FR_2.4.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_2.4.5_rarefaction_curves.qzv \
--p-min-depth 900 \
--p-max-depth 152786

  ## Core-metrics
qiime diversity core-metrics \
--i-table sidle/reconstruction/FR_2.4.5_tabla_recon.qza \
--m-metadata-file metadata_mama.tsv \
--p-sampling-depth 20054 \
--output-dir sidle/reconstruction/diversidad/FR_2.4.5

  ## Diversidad Alfa
qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_2.4.5/observed_features_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_2.4.5/observed_features_vector.qzv

qiime diversity alpha-group-significance \
--i-alpha-diversity sidle/reconstruction/diversidad/FR_2.4.5/shannon_vector.qza \
--m-metadata-file metadata_mama.tsv \
--o-visualization sidle/reconstruction/diversidad/FR_2.4.5/shannon_vector.qzv

  ## Diversidad Beta
qiime diversity beta-group-significance \
--i-distance-matrix sidle/reconstruction/diversidad/FR_2.4.5/bray_curtis_distance_matrix.qza \
--m-metadata-file metadata_mama.tsv \
--m-metadata-column Subtype \
--o-visualization sidle/reconstruction/diversidad/FR_2.4.5/bray_curtis_distance_matrix_Subtype.qzv \
--p-pairwise   
