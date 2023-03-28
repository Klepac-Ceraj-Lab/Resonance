ordered_brain_segments_list = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-cerebellum-white-matter", "right-cerebellum-white-matter",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-caudate", "right-caudate",
    "left-putamen", "right-putamen",
    "left-pallidum", "right-pallidum",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-basal-forebrain", "right-basal-forebrain", 
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lingual", "right-lingual",
    "left-parahippocampal", "right-parahippocampal",
    "left-paracentral", "right-paracentral",
    "left-pars-opercularis", "right-pars-opercularis",
    "left-pars-orbitalis", "right-pars-orbitalis",
    "left-pars-triangularis", "right-pars-triangularis",
    "left-pericalcarine", "right-pericalcarine",
    "left-postcentral", "right-postcentral",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-precentral", "right-precentral",
    "left-precuneus", "right-precuneus",
    "left-superior-frontal", "right-superior-frontal",
    "left-supramarginal", "right-supramarginal",
    "left-insula", "right-insula",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
    "Brain-stem", "CSF"
]

interesting_segments = [
    "left-temporal", "right-temporal",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-parietal", "right-parietal",
    "left-middle-frontal", "right-middle-frontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-hippocampus", "right-hippocampus",
    "left-amygdala", "right-amygdala",
    "left-accumbens-area", "right-accumbens-area",
    "left-cuneus", "right-cuneus", 
    "left-entorhinal", "right-entorhinal",
    "left-fusiform", "right-fusiform",
    "left-isthmus-cingulate", "right-isthmus-cingulate",
    "left-lingual", "right-lingual",
    "left-parahippocampal", "right-parahippocampal",
    "left-paracentral", "right-paracentral",
    "left-pars-opercularis", "right-pars-opercularis",
    "left-pars-orbitalis", "right-pars-orbitalis",
    "left-pericalcarine", "right-pericalcarine",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-precentral", "right-precentral",
    "left-precuneus", "right-precuneus",
    "left-superior-frontal", "right-superior-frontal",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
]

symmetric_segment_strings = [
    "temporal", "temporal",
    "orbitofrontal", "orbitofrontal",
    "parietal", "parietal",
    "middle-frontal", "middle-frontal",
    "anterior-cingulate", "anterior-cingulate",
    "lateral-occipital", "lateral-occipital",
    "thalamus-proper", "thalamus-proper",
    "hippocampus", "hippocampus",
    "amygdala", "amygdala",
    "accumbens-area", "accumbens-area",
    "cuneus", "cuneus", 
    "entorhinal", "entorhinal",
    "fusiform", "fusiform",
    "isthmus-cingulate", "isthmus-cingulate",
    "lingual", "lingual",
    "parahippocampal", "parahippocampal",
    "paracentral", "paracentral",
    "pars-opercularis", "pars-opercularis",
    "pars-orbitalis", "pars-orbitalis",
    "pericalcarine", "pericalcarine",
    "posterior-cingulate", "posterior-cingulate",
    "precentral", "precentral",
    "precuneus", "precuneus",
    "superior-frontal", "superior-frontal",
    "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
]

interesting_taxa = [
    "Adlercreutzia_equolifaciens",
    "Agathobaculum_butyriciproducens",
    "Akkermansia_muciniphila",
    "Alistipes_finegoldii",
    "Alistipes_putredinis",
    "Anaerostipes_hadrus",
    "Asaccharobacter_celatus",
    "Bacteroides_caccae",
    "Bacteroides_ovatus",
    "Bacteroides_stercoris",
    "Bacteroides_thetaiotaomicron",
    "Bacteroides_uniformis",
    "Bacteroides_vulgatus",
    "Bifidobacterium_longum",
    "Bifidobacterium_pseudocatenulatum",
    "Blautia_obeum",
    "Blautia_wexlerae",
    "Collinsella_aerofaciens",
    "Coprococcus_comes",
    "Coprococcus_eutactus",
    "Dorea_formicigenerans",
    "Dorea_longicatena",
    "Escherichia_coli",
    "Eubacterium_eligens",
    "Eubacterium_hallii",
    "Eubacterium_rectale",
    "Faecalibacterium_prausnitzii",
    "Flavonifractor_plautii",
    "Fusicatenibacter_saccharivorans",
    "Gemmiger_formicilis",
    "Intestinibacter_bartlettii",
    "Parabacteroides_distasonis",
    "Parasutterella_excrementihominis",
    "Prevotella_copri",
    "Roseburia_faecis",
    "Roseburia_intestinalis",
    "Roseburia_inulinivorans",
    "Ruminococcus_bicirculans",
    "Ruminococcus_gnavus",
    "Ruminococcus_torques",
    "Streptococcus_salivarius",
    "Streptococcus_thermophilus"
]