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


# interesting_segments = [ ## OLD LIST
#     "left-temporal", "right-temporal",
#     "left-orbitofrontal", "right-orbitofrontal",
#     "left-parietal", "right-parietal",
#     "left-middle-frontal", "right-middle-frontal",
#     "left-anterior-cingulate", "right-anterior-cingulate",
#     "left-lateral-occipital", "right-lateral-occipital",
#     "left-thalamus-proper", "right-thalamus-proper",
#     "left-hippocampus", "right-hippocampus",
#     "left-amygdala", "right-amygdala",
#     "left-accumbens-area", "right-accumbens-area",
#     "left-cuneus", "right-cuneus", 
#     "left-entorhinal", "right-entorhinal",
#     "left-fusiform", "right-fusiform",
#     "left-isthmus-cingulate", "right-isthmus-cingulate",
#     "left-lingual", "right-lingual",
#     "left-parahippocampal", "right-parahippocampal",
#     "left-paracentral", "right-paracentral",
#     "left-pars-opercularis", "right-pars-opercularis",
#     "left-pars-orbitalis", "right-pars-orbitalis",
#     "left-pericalcarine", "right-pericalcarine",
#     "left-posterior-cingulate", "right-posterior-cingulate",
#     "left-precentral", "right-precentral",
#     "left-precuneus", "right-precuneus",
#     "left-superior-frontal", "right-superior-frontal",
#     "cerebellar-vermal-lobules-I-V", "cerebellar-vermal-lobules-VI-VII", "cerebellar-vermal-lobules-VIII-X",
# ]

# From Figures of Merit:
# postcentral
# pars-triangularis
# entorhinal
# amygdala
# caudate
# precuneus
# thalamus-proper
# parahippocampal
# lateral-occipital
# orbitofrontal
# fusiform
# anterior-cingulate
# pallidum
# precentral
# paracentral
# pericalcarine
# accumbens-area
# cuneus
# cerebellum-white-matter - super loaded on Age
# lingual

interesting_segments = [
    "left-lingual", "right-lingual",
    "left-cuneus", "right-cuneus", 
    "left-accumbens-area", "right-accumbens-area",
    "left-pericalcarine", "right-pericalcarine",
    "left-precentral", "right-precentral",
    "left-paracentral", "right-paracentral",
    "left-pallidum", "right-pallidum",
    "left-fusiform", "right-fusiform",
    "left-amygdala", "right-amygdala",
    "left-orbitofrontal", "right-orbitofrontal",
    "left-anterior-cingulate", "right-anterior-cingulate",
    "left-posterior-cingulate", "right-posterior-cingulate",
    "left-lateral-occipital", "right-lateral-occipital",
    "left-thalamus-proper", "right-thalamus-proper",
    "left-entorhinal", "right-entorhinal",
    "left-parahippocampal", "right-parahippocampal",
    "left-precuneus", "right-precuneus"
]

symmetric_segment_strings =  replace.(interesting_segments, r"(left|right)-"=>"")

interesting_taxa = [
    # "Adlercreutzia_equolifaciens", # Not present anymore after Chocophlan version update
    # "Agathobaculum_butyriciproducens", # Removed on revision pruning
    "Akkermansia_muciniphila",
    "Alistipes_finegoldii",
    # "Alistipes_putredinis", # Removed on revision pruning
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
    # "Blautia_obeum", # Removed on revision pruning
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
    # "Gemmiger_formicilis", # Removed on revision pruning
    # "Intestinibacter_bartlettii", # Removed on revision pruning
    # "Parabacteroides_distasonis", # Removed on revision pruning
    "Parasutterella_excrementihominis",
    "Prevotella_copri",
    # "Roseburia_faecis", # Removed on revision pruning
    "Roseburia_intestinalis",
    # "Roseburia_inulinivorans", # Removed on revision pruning
    # "Ruminococcus_bicirculans", # Removed on revision pruning
    # "Ruminococcus_gnavus", # Removed on revision pruning
    "Ruminococcus_torques",
    "Streptococcus_salivarius",
    "Streptococcus_thermophilus"
]