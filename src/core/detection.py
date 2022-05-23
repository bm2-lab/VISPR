from .processor import Processor
import os
import json

class Detection:
    def __init__(self,
                 strg_gtf,
                 strg_gtf_fpkm_filter,
                 ref_fa_fixed,
                 viral_gtf,
                 viral_ref_fa,
                 chrom_name,
                 output_dir,
                 bl_xlsx,
                 name,
                 bowtie_nprocesses,
                 bowtie_nthreads,
                 bowtie_idx_prefix,
                 bowtie_bg_rna,
                 virus_gene_selected,
                 cas_definition_json,
                 grna_mismatch):
        self.strg_gtf = strg_gtf
        self.strg_gtf_fpkm_filter = strg_gtf_fpkm_filter
        self.ref_fa_fixed = ref_fa_fixed
        self.viral_gtf = viral_gtf
        self.viral_ref_fa = viral_ref_fa
        self.chrom_name = chrom_name
        self.output_dir = output_dir
        self.bl_xlsx = bl_xlsx
        self.name = name
        self.bowtie_nprocesses = bowtie_nprocesses
        self.bowtie_nthreads = bowtie_nthreads
        self.virus_gene_selected = virus_gene_selected
        if cas_definition_json is None:
            self.cas_definition = Processor.g_default_cas_definition
        else:
            self.cas_definition = json.load(open(cas_definition_json, encoding='utf-8'))
        self.grna_mismatch = grna_mismatch

        os.makedirs(self.output_dir, exist_ok=True)
        self.viral_gtf_json = f'{self.output_dir}/{os.path.splitext(os.path.basename(self.viral_gtf))[0]}.json'
        self.bowtie_out_dir = f'{self.output_dir}/bowtie_out'
        os.makedirs(self.bowtie_out_dir, exist_ok=True)
        
        if bowtie_idx_prefix is None:
            if not os.path.exists(self.strg_gtf):
                raise RuntimeError(f'GTF file {self.strg_gtf} not exist')
            if not os.path.exists(self.ref_fa_fixed):
                raise RuntimeError(f'Reference file {self.ref_fa_fixed} not exist')
            self.strg_gtf_json = f'{self.output_dir}/{os.path.splitext(os.path.basename(self.strg_gtf))[0]}.json'
            
            self.bg_rna_fa = f'{self.output_dir}/{os.path.splitext(os.path.basename(self.ref_fa_fixed))[0]}_bg_rna.fa'
            self.bg_rna_idx_dir = f'{self.bowtie_out_dir}/bowtie2_idx'
            os.makedirs(self.bg_rna_idx_dir, exist_ok=True)
            self.bg_rna_idx_prefix = f'{self.bg_rna_idx_dir}/{self.name}_bg_rna'
            self.idx_presetted = False
        else:
            self.bg_rna_fa = bowtie_bg_rna
            self.bg_rna_idx_prefix = bowtie_idx_prefix
            self.idx_presetted = True

        self.cas_dt = {cas:(self.cas_definition[cas]['pfs'], self.cas_definition[cas]['loc']) for cas in self.cas_definition}
        self.bowtie_prefix = {}
        for cas in self.cas_definition:
            for l in range(self.cas_definition[cas]['min_len'], self.cas_definition[cas]['max_len']+1):
                os.makedirs(f'{self.bowtie_out_dir}/{cas}/{self.name}_{cas}_{l}', exist_ok=True)
                self.bowtie_prefix[(cas, l)] = f'{self.bowtie_out_dir}/{cas}/{self.name}_{cas}_{l}/{self.name}_{cas}_{l}'

    def run(self):
        if not self.idx_presetted:
            Processor.analyze_strg_gtf(self.strg_gtf, self.ref_fa_fixed, self.strg_gtf_json, self.strg_gtf_fpkm_filter, self.name)

            Processor.build_background_rna_index(self.strg_gtf_json, self.bg_rna_fa, self.bg_rna_idx_prefix, self.bowtie_nthreads)

        Processor.analyze_viral_gtf(self.viral_gtf, self.viral_ref_fa, self.viral_gtf_json, self.name)

        Processor.generate_virus_grna(self.viral_gtf_json, self.virus_gene_selected, self.bowtie_prefix, self.cas_dt, self.chrom_name)

        Processor.mapping_virus_grna(self.bowtie_prefix, self.bowtie_out_dir, self.bg_rna_idx_prefix, self.grna_mismatch, self.bowtie_nthreads, self.bowtie_nprocesses)
        
        if self.bl_xlsx:
            Processor.detection_report(self.bowtie_prefix, self.cas_dt, self.output_dir, self.name, self.bg_rna_fa)

    