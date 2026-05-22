# ARI_patients-metagenomics
Closing the diagnostic gap in respiratory surveillance: Metagenomic characterization of pathogens in undiagnosed ARI patients

# 1. DATABASES CONFIGURATION AND OUTPUT DIRECTORIES
```r
BOWTIE2_INDEX="/data/iris/GRCh38_noalt_as/GRCh38_noalt_as"
HISAT2_INDEX="/data/iris/grch38/genome"
KRAKEN2_DB="/data/iris/kraken_db/viral"

mkdir -p cleaned mapped unmapped assemblies kraken kraken_r krona
```

# 2. GENERAL LOOP
```r
for R1 in *_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="${SAMPLE}_R2_001.fastq.gz"

    echo "? Processando $SAMPLE..."

    ## 1. fastp - cleaning
    echo "running fastp"
    fastp -i $R1 -I $R2 \
          -o cleaned/${SAMPLE}_clean_R1.fastq.gz \
          -O cleaned/${SAMPLE}_clean_R2.fastq.gz \
	  --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
          --qualified_quality_phred 20 \
          --detect_adapter_for_pe \
	  --correction  \
	  -h cleaned/${SAMPLE}_fastp.html \
          -j cleaned/${SAMPLE}_fastp.json

    ## 2.  bowtie2 - Removing human contamination
    echo "running bowtie2"
    bowtie2 -x $BOWTIE2_INDEX \
        -1 cleaned/${SAMPLE}_clean_R1.fastq.gz \
        -2 cleaned/${SAMPLE}_clean_R2.fastq.gz \
        --very-sensitive -p 8 \
        --un-conc-gz unmapped/${SAMPLE}_bowtie2_unmapped.fastq.gz \
        -S mapped/${SAMPLE}_bowtie2_mapped.sam

    ## 2.1. hisat2 (alternative)
     echo "running hisat2"
     hisat2 -x $HISAT2_INDEX \
         -1 unmapped/${SAMPLE}_bowtie2_unmapped.fastq.1.gz \
         -2 unmapped/${SAMPLE}_bowtie2_unmapped.fastq.2.gz \
         --un-conc-gz unmapped/${SAMPLE}_hisat2_unmapped.fastq.gz \
         -S mapped/${SAMPLE}_hisat2_mapped.sam
	 
     echo "rename files"
	 mv unmapped/${SAMPLE}_hisat2_unmapped.fastq.1.gz unmapped/${SAMPLE}_hisat2_unmapped.1.fastq.gz
         mv unmapped/${SAMPLE}_hisat2_unmapped.fastq.2.gz unmapped/${SAMPLE}_hisat2_unmapped.2.fastq.gz

   ## 3. kraken2 - reads taxonomic identification
    echo "running kraken2_r"
    kraken2 --db $KRAKEN2_DB \
        --paired unmapped/${SAMPLE}_hisat2_unmapped.1.fastq.gz \
                unmapped/${SAMPLE}_hisat2_unmapped.2.fastq.gz \
       --report kraken_r/${SAMPLE}_report.txt \
        --output kraken_r/${SAMPLE}_output.txt \
        --use-names
	 
	 
    ## 4. spades - assembly 
     
        echo "running MetaSPAdes"
    /home/iris/SPAdes-4.2.0-Linux/bin/metaspades.py -1 unmapped/${SAMPLE}_hisat2_unmapped.1.fastq.gz \
              -2 unmapped/${SAMPLE}_hisat2_unmapped.2.fastq.gz \
              -o assemblies/${SAMPLE}_spades \
              -t 8 
	      
    ## 4.1 kraken2 - taxonomic identification
    echo "running kraken2"
    kraken2 --db $KRAKEN2_DB \
        --report kraken/${SAMPLE}_report.txt \
        --output kraken/${SAMPLE}_output.txt \
        --use-names \
	assemblies/${SAMPLE}_spades/contigs.fasta 
	
	      
    ## 5. krona - visualization
    #cut -f2,3 kraken/${SAMPLE}_output.txt | \
        #ktImportTaxonomy -o krona/${SAMPLE}_krona.html -

    echo "? $SAMPLE finalizado."
done

echo "? Pipeline completo!"
```

# 3. KRAKENTOOLS
```r
# 3.1. environment activation 
conda activate krakentools

# 3.2. based on samples with  "_report.txt" termination
for report in *_report.txt; do
    
    # Extract prefix 
    SAMPLE=$(basename "$report" _report.txt)
    
    echo "Procesando muestra: $SAMPLE"
    
    # Assign output Kraken2 file names
    OUTPUT_KRAKEN="${SAMPLE}_output.txt"
    
    # Input name files configuration
    
    R1="${SAMPLE}_R1.fastq.gz"
    R2="${SAMPLE}_R2.fastq.gz"
    
    # Assign output FASTQ file names
    OUT_R1="lecturas_${SAMPLE}_R1.fastq"
    OUT_R2="lecturas_${SAMPLE}_R2.fastq"
    
    # run only if input files are present
    if [[ -f "$R1" && -f "$R2" ]]; then
        extract_kraken_reads.py \
          -k "$OUTPUT_KRAKEN" \
          -s1 "$R1" \
          -s2 "$R2" \
          -t 297248 \
          -o "$OUT_R1" \
          -o2 "$OUT_R2" \
          -r "$report" \
          --include-children \
          --fastq-output
    else
        echo "Error: No se encontraron los archivos FASTQ para $SAMPLE"
    fi

done

echo "all samples has been proccessed"

```

# 4. BOWTIE2
```r

# 1. Generate indexes - BOWTIE2

if [ ! -f "idx.1.bt2" ]; then
    bowtie2-build ref.fasta idx
else
    echo "El índice idx ya existe. Saltando paso..."
fi

# Nombre del reporte final de cobertura
REPORT_FINAL="resumen_cobertura.txt"
echo "REPORTE DE COBERTURA Y PROFUNDIDAD enterovirusA ===" > "$REPORT_FINAL"


# 2. LOOP
for R1 in lecturas_*_R1.fastq; do
    
    # Verificar si existen archivos
    [ -e "$R1" ] || { echo " No se encontraron archivos de lectura R1 con el patrón 'lecturas_*_R1.fastq'."; exit 1; }

    # Extraer el nombre de la muestra (borra 'lecturas_' y '_R1.fastq')
    SAMPLE=$(basename "$R1" _R1.fastq | sed 's/lecturas_//')
    
    # Definir el archivo R2 correspondiente
    R2="lecturas_${SAMPLE}_R2.fastq"
    
    echo " PROCESANDO MUESTRA: $SAMPLE"

    if [ ! -f "$R2" ]; then
        echo " Error: No se encontró el par R2 ($R2). Saltando..."
        continue
    fi

    # Mapeo con Bowtie2 usando la configuración muy sensible
    echo "[1/4] Ejecutando Bowtie2 (--very-sensitive)..."
    bowtie2 --very-sensitive -x idx -1 "$R1" -2 "$R2" -S "${SAMPLE}.sam"

    # Conversión y ordenamiento de SAM a BAM
    echo "[2/4] Ordenando y convirtiendo SAM a BAM..."
    samtools sort -o "${SAMPLE}.sorted.bam" "${SAMPLE}.sam"

    # Creación del índice del archivo BAM (.bai)
    echo "[3/4] Creando índice .bai..."
    samtools index "${SAMPLE}.sorted.bam"

    # Cálculo de la cobertura horizontal y profundidad
    echo "[4/4] Calculando cobertura..."
    samtools coverage "${SAMPLE}.sorted.bam" >> "$REPORT_FINAL"

    # Limpieza del archivo SAM intermedio para ahorrar espacio
    rm "${SAMPLE}.sam"
done
echo "DONE"
```

# 5. RPM 
```r
# Taxon ID oficial de respirovirus
TAXON_ID="XXXX" 

# Archivo de salida estructurado (formato TSV compatible con Excel)
OUTPUT_TABLA="resultados_carga_viral.txt"

# Crear el encabezado del reporte
echo -e "Muestra\tLecturas_Totales_Nt\tLecturas\tRPM_Teorica" > "$OUTPUT_TABLA"

# El bucle busca todos tus archivos _output.txt de Kraken en la carpeta
for kraken_out in *_output.txt; do
    
    # Comprobar la existencia de archivos en la carpeta
    [ -e "$kraken_out" ] || { echo "No se encontraron archivos con terminación _output.txt"; exit 1; }
    
    # Extraer el nombre de la muestra limpiando el sufijo
    SAMPLE=$(basename "$kraken_out" _output.txt)
    
    # 1. Número total de lecturas secuenciadas (Nt)
    NT=$(wc -l < "$kraken_out")
    
    # Validar si el archivo está vacío para prevenir un error de división por cero
    if [ "$NT" -eq 0 ]; then
        echo " Muestra $SAMPLE vacía (Nt = 0). Saltando..."
        continue
    fi
    
    # 2. Número de lecturas asignadas a enterovirusA (Nv)
    # Buscamos el Taxon ID en la tercera columna o mapeos taxonómicos secundarios del formato Kraken
    NV=$(awk -v tax="$TAXON_ID" '$3 == tax || $0 ~ tax' "$kraken_out" | wc -l)
    
    # 3. Cálculo automatizado de la Carga Viral con tu comando de Python
    RPM=$(python3 -c "print(f'{($NV * 1000000) / $NT:.2f}')")
    
    # Mostrar resultados en tiempo real en la pantalla de la terminal
    echo "Muestra: $SAMPLE | Nt: $NT | Nv: $NV | RPM: $RPM"
    
    # Guardar los datos en la tabla acumulativa
    echo -e "$SAMPLE\t$NT\t$NV\t$RPM" >> "$OUTPUT_TABLA"

done
```

# 6. SMRN
```r
# Archivo de salida final
OUTPUT_SMRN="resultados_smrn.txt"

# Crear encabezado de la tabla
echo -e "Muestra\tSMRN_Lecturas_Estrictas" > "$OUTPUT_SMRN"

# El bucle busca todos los archivos BAM ordenados
for bam_file in *.sorted.bam; do
    
    # Comprobar si existen archivos
    [ -e "$bam_file" ] || { echo "No se encontraron archivos .sorted.bam"; exit 1; }
    
    # Extraer el nombre de la muestra (ej: de .sorted.bam toma S2)
    SAMPLE=$(basename "$bam_file" .sorted.bam)
    
    # Calcular el SMRN usando tu configuración de Samtools
    # Nota: Funciona igual en archivos .bam (que son más rápidos de procesar)
    SMRN=$(samtools view -c -f 1 -F 12 "$bam_file")
    
    # Mostrar el progreso en la terminal
    echo " Muestra: $SAMPLE | SMRN: $SMRN"
    
    # Guardar en el archivo de texto permanente
    echo -e "$SAMPLE\t$SMRN" >> "$OUTPUT_SMRN"

done

echo "DONE"
```
