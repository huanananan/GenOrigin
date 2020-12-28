import os
import json


def read_json(homologyB_dir):
    file = open(homologyB_dir, 'r', encoding='UTF-8')
    homologyB = json.loads(file.read())
    file.close()
    return homologyB


work_dir = 'GenOrigin/tree/master/'
division = 'vertebrates'
species = 'homo_sapiens'
homologyB = read_json(work_dir + "homology_json/" + division + "/" + species)
for i in range(int(len(homologyB[species]) / 500) + 1):
    python_file1 = """import os
import subprocess
import json
from Bio import Phylo
from collections import defaultdict

analyse_group = {
    "metazoa": [["Drosophila_mojavensis", "Octopus_bimaculoides", "Tribolium_castaneum", "Atta_cephalotes",
                 "Acyrthosiphon_pisum",
                 "Caenorhabditis_briggsae", "Drosophila_ananassae", "Rhodnius_prolixus", "Brugia_malayi",
                 "Drosophila_yakuba",
                 "Stegodyphus_mimosarum", "Heliconius_melpomene", "Caenorhabditis_elegans", "Anopheles_darlingi",
                 "Zootermopsis_nevadensis", "Drosophila_simulans", "Orchesella_cincta", "Culex_quinquefasciatus",
                 "Mayetiola_destructor", "Solenopsis_invicta", "Melitaea_cinxia", "Belgica_antarctica",
                 "Caenorhabditis_japonica",
                 "Caenorhabditis_remanei", "Anopheles_gambiae", "Lepeophtheirus_salmonis", "Crassostrea_gigas",
                 "Dendroctonus_ponderosae", "Branchiostoma_lanceolatum", "Pediculus_humanus", "Pristionchus_pacificus",
                 "Apis_mellifera", "Strigamia_maritima", "Amphimedon_queenslandica", "Mnemiopsis_leidyi",
                 "Drosophila_erecta",
                 "Drosophila_sechellia", "Danaus_plexippus", "Trichoplax_adhaerens", "Teleopsis_dalmanni",
                 "Daphnia_pulex",
                 "Bombyx_mori", "Daphnia_magna", "Aedes_aegypti", "Drosophila_pseudoobscura",
                 "Leptotrombidium_deliense",
                 "Drosophila_grimshawi", "Drosophila_persimilis", "Helobdella_robusta", "Lucilia_cuprina",
                 "Drosophila_melanogaster", "Drosophila_willistoni", "Caenorhabditis_brenneri", "Drosophila_virilis",
                 "Nematostella_vectensis", "Capitella_teleta"]],

    "fungi": [["Eutypa_lata", "Armillaria_ostoyae", "Conidiobolus_coronatus", "Trichoderma_citrinoviride",
               "Torulaspora_delbrueckii", "Byssochlamys_spectabilis", "Pisolithus_tinctorius",
               "Fusarium_fujikuroi", "Ophiostoma_piceae", "Coniophora_puteana", "Debaryomyces_hansenii",
               "Diplodia_seriata", "Phaeosphaeria_nodorum", "Pseudocercospora_fijiensis",
               "Tetrapisispora_blattae", "Cordyceps_militaris", "Exophiala_dermatitidis",
               "Jaapia_argillacea", "Macrophomina_phaseolina", "Colletotrichum_graminicola",
               "Plicaturopsis_crispa", "Punctularia_strigosozonata", "Blumeria_graminis",
               "Clavispora_lusitaniae", "Saccharomyces_cerevisiae", "Aspergillus_fumigatus",
               "Alternaria_alternata", "Botrytis_cinerea", "Thermothelomyces_thermophila",
               "Aspergillus_flavus", "Trichoderma_longibrachiatum", "Corynespora_cassiicola",
               "Anthracocystis_flocculosa", "Scheffersomyces_stipitis", "Aspergillus_clavatus",
               "Pleurotus_ostreatus", "Aspergillus_terreus", "Sporisorium_reilianum",
               "Schizosaccharomyces_pombe", "Blastomyces_dermatitidis", "Leptosphaeria_maculans",
               "Coccidioides_immitis", "Batrachochytrium_dendrobatidis", "Heterobasidion_irregulare",
               "Tuber_melanosporum", "Trichoderma_harzianum", "Fusarium_verticillioides",
               "Colletotrichum_gloeosporioides", "Neurospora_tetrasperma", "Drechslerella_stenobrocha",
               "Spizellomyces_punctatus", "Aspergillus_oryzae", "Gelatoporia_subvermispora",
               "Phaeoacremonium_minimum", "Valsa_mali", "Dactylellina_haptotyla",
               "Zygosaccharomyces_rouxii", "Sphaerulina_musiva", "Zymoseptoria_tritici",
               "Naumovozyma_dairenensis", "Coniochaeta_ligniaria", "Claviceps_purpurea",
               "Thielavia_terrestris", "Trichophyton_benhamiae", "Agaricus_bisporus", "Stereum_hirsutum",
               "Tuber_borchii", "Fusarium_oxysporum", "Ustilago_hordei", "Paracoccidioides_brasiliensis",
               "Ascosphaera_apis", "Kazachstania_naganishii", "Trichoderma_atroviride",
               "Candida_tropicalis", "Trichoderma_asperellum", "Tolypocladium_paradoxum",
               "Postia_placenta", "Verticillium_dahliae", "Neurospora_crassa", "Erysiphe_necator",
               "Pneumocystis_carinii", "Tetrapisispora_phaffii", "Serendipita_vermifera",
               "Bipolaris_maydis", "Brettanomyces_bruxellensis", "Allomyces_macrogynus", "Tuber_magnatum",
               "Stachybotrys_chartarum", "Rhizophagus_clarus", "Tolypocladium_ophioglossoides",
               "Aspergillus_niger", "Armillaria_gallica", "Laccaria_bicolor", "Calocera_cornea",
               "Kazachstania_africana", "Neolentinus_lepideus", "Colletotrichum_higginsianum",
               "Serendipita_indica", "Ustilago_maydis", "Diplodia_corticola", "Metarhizium_anisopliae",
               "Fomitopsis_pinicola", "Hydnomerulius_pinastri", "Fomitiporia_mediterranea",
               "Coniosporium_apollinis", "Ganoderma_sinense", "Metarhizium_album", "Madurella_mycetomatis",
               "Wolfiporia_cocos", "Phanerochaete_carnosa", "Metarhizium_rileyi", "Candida_albicans",
               "Kluyveromyces_lactis", "Trichoderma_virens", "Serpula_lacrymans",
               "Vanderwaltozyma_polyspora", "Meyerozyma_guilliermondii", "Aspergillus_nidulans",
               "Cryptococcus_neoformans", "Histoplasma_capsulatum", "Komagataella_pastoris",
               "Naumovozyma_castellii", "Tremella_mesenterica", "Purpureocillium_lilacinum",
               "Tolypocladium_capitatum", "Metarhizium_acridum", "Pochonia_chlamydosporia",
               "Coccidioides_posadasii", "Eremothecium_gossypii", "Pyrenophora_triticirepentis",
               "Trichoderma_reesei", "Scleroderma_citrinum", "Cordyceps_confragosa",
               "Trichophyton_rubrum"]],

    "vertebrates": [["Cyprinodon_variegatus", "Jaculus_jaculus", "Callithrix_jacchus", "Panthera_pardus",
                     "Cavia_aperea", "Homo_sapiens", "Aotus_nancymaae", "Scophthalmus_maximus",
                     "Mus_spicilegus", "Meriones_unguiculatus", "Gallus_gallus", "Takifugu_rubripes",
                     "Meleagris_gallopavo", "Poecilia_reticulata", "Manacus_vitellinus", "Pan_paniscus",
                     "Bos_taurus", "Gadus_morhua", "Anolis_carolinensis", "Cebus_capucinus",
                     "Rattus_norvegicus", "Dipodomys_ordii", "Ovis_aries", "Vulpes_vulpes",
                     "Mandrillus_leucophaeus", "Ficedula_albicollis", "Propithecus_coquereli",
                     "Otolemur_garnettii", "Ursus_maritimus", "Lonchura_striata",
                     "Dromaius_novaehollandiae", "Gasterosteus_aculeatus", "Vombatus_ursinus",
                     "Lates_calcarifer", "Pelodiscus_sinensis", "Monopterus_albus",
                     "Electrophorus_electricus", "Dasypus_novemcinctus", "Betta_splendens",
                     "Mustela_putorius", "Junco_hyemalis", "Lepidothrix_coronata",
                     "Erpetoichthys_calabaricus", "Anser_brachyrhynchus", "Macaca_nemestrina",
                     "Ictidomys_tridecemlineatus", "Ciona_intestinalis", "Esox_lucius", "Chrysemys_picta",
                     "Cavia_porcellus", "Vicugna_pacos", "Calidris_pygmaea", "Hippocampus_comes",
                     "Sorex_araneus", "Equus_asinus", "Zonotrichia_albicollis", "Sarcophilus_harrisii",
                     "Carlito_syrichta", "Cyanistes_caeruleus", "Pan_troglodytes", "Mola_mola",
                     "Spermophilus_dauricus", "Pygocentrus_nattereri", "Saimiri_boliviensis",
                     "Urocitellus_parryii", "Rhinopithecus_roxellana", "Mus_spretus",
                     "Nomascus_leucogenys", "Melopsittacus_undulatus", "Amphiprion_ocellaris",
                     "Macaca_mulatta", "Amphiprion_percula", "Taeniopygia_guttata", "Colobus_angolensis",
                     "Bos_mutus", "Maylandia_zebra", "Notamacropus_eugenii", "Stegastes_partitus",
                     "Eptatretus_burgeri", "Anabas_testudineus", "Cottoperca_gobio",
                     "Ornithorhynchus_anatinus", "Cricetulus_griseus", "Mesocricetus_auratus",
                     "Xiphophorus_maculatus", "Echinops_telfairi", "Marmota_marmota", "Myotis_lucifugus",
                     "Loxodonta_africana", "Labrus_bergylta", "Oryzias_latipes", "Felis_catus",
                     "Gopherus_agassizii", "Coturnix_japonica", "Sus_scrofa", "Panthera_tigris",
                     "Monodelphis_domestica", "Pteropus_vampyrus", "Denticeps_clupeoides", "Canis_lupus",
                     "Apteryx_owenii", "Pogona_vitticeps", "Apteryx_rowi", "Notechis_scutatus",
                     "Ictalurus_punctatus", "Piliocolobus_tephrosceles", "Numida_meleagris",
                     "Neovison_vison", "Peromyscus_maniculatus", "Parambassis_ranga", "Poecilia_mexicana",
                     "Parus_major", "Bos_indicus", "Neolamprologus_brichardi", "Phascolarctos_cinereus",
                     "Ochotona_princeps", "Latimeria_chalumnae", "Gambusia_affinis", "Sphenodon_punctatus",
                     "Larimichthys_crocea", "Erinaceus_europaeus", "Macaca_fascicularis",
                     "Salvator_merianae", "Fundulus_heteroclitus", "Cercocebus_atys",
                     "Chlorocebus_sabaeus", "Pundamilia_nyererei", "Papio_anubis", "Equus_caballus",
                     "Ailuropoda_melanoleuca", "Seriola_dumerili", "Xiphophorus_couchianus", "Bison_bison",
                     "Hucho_hucho", "Fukomys_damarensis", "Crocodylus_porosus", "Tetraodon_nigroviridis",
                     "Amphilophus_citrinellus", "Microtus_ochrogaster", "Nannospalax_galili",
                     "Anas_platyrhynchos", "Microcebus_murinus", "Nothoprocta_perdicaria",
                     "Gorilla_gorilla", "Mus_musculus", "Choloepus_hoffmanni", "Chelonoidis_abingdonii",
                     "Callorhinchus_milii", "Octodon_degus", "Scleropages_formosus", "Ursus_americanus",
                     "Mastacembelus_armatus", "Rhinopithecus_bieti", "Haplochromis_burtoni", "Mus_caroli",
                     "Gouania_willdenowi", "Serinus_canaria", "Castor_canadensis", "Heterocephalus_glaber",
                     "Procavia_capensis", "Petromyzon_marinus", "Clupea_harengus",
                     "Acanthochromis_polyacanthus", "Tupaia_belangeri", "Capra_hircus",
                     "Theropithecus_gelada", "Chinchilla_lanigera", "Lepisosteus_oculatus",
                     "Calidris_pugnax", "Apteryx_haastii", "Prolemur_simus", "Xenopus_tropicalis",
                     "Oryctolagus_cuniculus", "Tursiops_truncatus", "Mus_pahari", "Danio_rerio",
                     "Pongo_abelii", "Seriola_lalandi"]],

    "protists": [
        ["Plasmodium_reichenowi", "Plasmodium_berghei", "Emiliania_huxleyi", "Phytophthora_infestans",
         "Phytophthora_sojae", "Thalassiosira_oceanica", "Aureococcus_anophagefferens",
         "Phytophthora_parasitica", "Chroomonas_mesostigmatica", "Phaeodactylum_tricornutum",
         "Phytophthora_ramorum", "Thalassiosira_pseudonana", "Ectocarpus_siliculosus", "Giardia_lamblia",
         "Capsaspora_owczarzaki", "Paramecium_tetraurelia", "Entamoeba_histolytica", "Plasmodium_chabaudi",
         "Eimeria_tenella", "Symbiodinium_microadriaticum", "Plasmodium_knowlesi", "Theileria_parva",
         "Perkinsus_marinus", "Sphaeroforma_arctica", "Plasmodium_inui", "Plasmodium_malariae",
         "Plasmodium_gonderi", "Dictyostelium_discoideum", "Leishmania_infantum", "Leishmania_major",
         "Thraustotheca_clavata", "Plasmodium_ovale", "Leishmania_donovani", "Naegleria_gruberi",
         "Reticulomyxa_filosa", "Plasmodiophora_brassicae", "Plasmodium_fragile", "Trypanosoma_brucei",
         "Plasmodium_falciparum", "Plasmodium_cynomolgi", "Toxoplasma_gondii", "Bigelowiella_natans",
         "Cryptosporidium_parvum", "Stylonychia_lemnae", "Hyaloperonospora_arabidopsidis",
         "Plasmodium_yoelii", "Albugo_laibachii", "Pseudonitzschia_multistriata",
         "Tetrahymena_thermophila", "Plasmodium_vivax", "Fragilariopsis_cylindrus", "Pythium_ultimum",
         "Saprolegnia_parasitica", "Plasmodium_coatneyi", "Plasmodium_gaboni", "Plasmodium_gallinaceum",
         "Nannochloropsis_gaditana"]],

    "bacteria": [
        ["Xanthomonas_campestris", "Lactococcus_lactis", "Buchnera_aphidicola", "Staphylococcus_aureus",
         "Geobacter_sulfurreducens", "Chlorobium_tepidum", "Fusobacterium_nucleatum",
         "Mannheimia_haemolytica", "Ralstonia_solanacearum", "Aggregatibacter_actinomycetemcomitans",
         "Mycobacterium_tuberculosis", "Gardnerella_vaginalis", "Bifidobacterium_longum",
         "Mesoplasma_florum", "Acinetobacter_baumannii", "Leuconostoc_mesenteroides",
         "Bordetella_pertussis", "Sinorhizobium_meliloti", "Flavobacterium_psychrophilum",
         "Shigella_dysenteriae", "Thermosynechococcus_elongatus", "Campylobacter_jejuni",
         "Moraxella_catarrhalis", "Microcystis_aeruginosa", "Coxiella_burnetii", "Haemophilus_influenzae",
         "Azotobacter_vinelandii", "Rhodopirellula_baltica", "Clostridioides_difficile",
         "Nostoc_punctiforme", "Thermus_thermophilus", "Bartonella_henselae", "Anaplasma_phagocytophilum",
         "Burkholderia_pseudomallei", "Bacteroides_thetaiotaomicron", "Propionibacterium_acnes",
         "Rhodobacter_sphaeroides", "Escherichia_coli", "Francisella_tularensis", "Legionella_pneumophila",
         "Micrococcus_luteus", "Myxococcus_xanthus", "Prevotella_intermedia", "Citrobacter_freundii",
         "Moorella_thermoacetica", "Rhodospirillum_rubrum", "Pseudomonas_aeruginosa", "Yersinia_pestis",
         "Helicobacter_pylori", "Lysinibacillus_sphaericus", "Shewanella_oneidensis",
         "Rickettsia_prowazekii", "Salmonella_enterica", "Paracoccus_denitrificans",
         "Listeria_monocytogenes", "Thermotoga_maritima", "Pasteurella_multocida",
         "Stenotrophomonas_maltophilia", "Gloeobacter_violaceus", "Actinobacillus_pleuropneumoniae",
         "Prochlorococcus_marinus", "Vibrio_cholerae", "Desulfovibrio_vulgaris", "Proteus_mirabilis",
         "Bacillus_subtilis", "Chloroflexus_aurantiacus", "Ureaplasma_parvum", "Mycoplasma_pneumoniae",
         "Porphyromonas_gingivalis", "Leptospira_interrogans", "Salinibacter_ruber",
         "Caulobacter_crescentus", "Enterococcus_faecalis", "Deinococcus_radiodurans",
         "Streptococcus_pneumoniae", "Vibrio_fischeri", "Corynebacterium_glutamicum",
         "Borrelia_burgdorferi", "Clostridium_botulinum", "Lactobacillus_plantarum", "Brucella_abortus",
         "Chlamydia_trachomatis"]],

    "plants": [["Chondrus_crispus", "Galdieria_sulphuraria", "Cyanidioschyzon_merolae"],
               ["Brassica_rapa", "Helianthus_annuus", "Brassica_oleracea", "Oryza_nivara", "Oryza_sativa",
                "Trifolium_pratense", "Eragrostis_tef", "Arabidopsis_thaliana", "Vitis_vinifera",
                "Coffea_canephora", "Prunus_persica", "Cucumis_sativus", "Leersia_perrieri",
                "Phaseolus_vulgaris", "Sorghum_bicolor", "Brassica_napus", "Capsicum_annuum",
                "Actinidia_chinensis", "Zea_mays", "Solanum_lycopersicum", "Lupinus_angustifolius",
                "Amborella_trichopoda", "Physcomitrella_patens", "Brachypodium_distachyon",
                "Oryza_brachyantha", "Triticum_aestivum", "Oryza_glumipatula", "Gossypium_raimondii",
                "Theobroma_cacao", "Aegilops_tauschii", "Triticum_urartu", "Musa_acuminata",
                "Arabidopsis_lyrata", "Hordeum_vulgare", "Arabidopsis_halleri", "Solanum_tuberosum",
                "Chlamydomonas_reinhardtii", "Oryza_meridionalis", "Oryza_glaberrima", "Vigna_radiata",
                "Triticum_turgidum", "Marchantia_polymorpha", "Populus_trichocarpa",
                "Ostreococcus_lucimarinus", "Oryza_barthii", "Oryza_longistaminata",
                "Selaginella_moellendorffii", "Daucus_carota", "Beta_vulgaris", "Glycine_max",
                "Medicago_truncatula", "Oryza_rufipogon", "Triticum_dicoccoides", "Oryza_punctata",
                "Manihot_esculenta"]]}

out_group = {"Ustilago_maydis": ["bacteria"], "Saccharomyces_cerevisiae": ["plants", "vertebrates"],
             "Homo_sapiens": ["plants"], "Drosophila_melanogaster": ["plants", "vertebrates"],
             "Caenorhabditis_elegans": ["plants", "vertebrates"]}


def read_json(homologyB_dir):
    file = open(homologyB_dir, 'r', encoding='UTF-8')
    homologyB = json.loads(file.read())
    file.close()
    return homologyB


def total_tree(gene_id, species, division, nwk_dir, trace, branch_length_dict, maximum_age, all_trace_dic, gene_age,
               gene_interval, gene_branch, ori_branch, gene_homology, homology_if_no_changed, species2pan_homology,
               pan_homology):
    try:
        species2pan_homology["_".join(species.split("_")[0:2]).capitalize()][gene_id]
    except:
        print('err 1')
        homology_changed = homology_if_no_changed
        return gene_age, gene_interval, gene_branch, ori_branch, gene_homology, homology_changed

    pan_gene_id = list(list(species2pan_homology[
                                "_".join(species.split("_")[0:2]).capitalize()][gene_id].values())[0].keys())[0]

    pan_species = list(list(species2pan_homology[
                                "_".join(species.split("_")[0:2]).capitalize()][gene_id].values())[0].values())[0][
        0].capitalize()

    try:
        pan_homology[species2species[pan_species]][pan_gene_id]
    except:
        print('err 2')
        homology_changed = homology_if_no_changed
        return gene_age, gene_interval, gene_branch, ori_branch, gene_homology, homology_changed

    for homology_species in pan_homology[species2species[pan_species]][pan_gene_id]:
        homology_if_no_changed[homology_species] = list(
            set(homology_if_no_changed[homology_species]) |
            set(pan_homology[species2species[pan_species]][pan_gene_id][homology_species])
        )
    homology_changed = homology_if_no_changed
    family2node2presence = run_count({species: {gene_id: homology_changed}}, nwk_dir, "-gain 1.4", species, division)

    gene_age, gene_interval, gene_branch, ori_branch, gene_homology = main_total(family2node2presence,
                                                                                 species, trace, branch_length_dict,
                                                                                 maximum_age, all_trace_dic, gene_id)

    return gene_age, gene_interval, gene_branch, ori_branch, gene_homology, homology_changed


def main_total(family2node2presence, species, trace, branch_length_dict, maximum_age, all_trace_dic,
               gene_id):
    # function

    tree = Phylo.read(nwk_dir, "newick")

    ori_branch = "_".join(species.split("_")[0:2]).capitalize()

    effect = []
    for node in family2node2presence[gene_id]:
        if family2node2presence[gene_id][node] > 0:
            if node.isdigit():
                effect.append('confidence=' + str(node))
            else:
                effect.append(node)

    trace_main_branch = {}
    for step in trace:
        if step in effect:
            ori_branch = step
            trace_main_branch[step] = 1
        else:
            break

    gene_homology = []
    for step in trace:
        gene_homology.append(step)
        if step not in trace_main_branch:
            trace_main_branch[step] = 0

    for species_homology in total_species_list:
        if family2node2presence[gene_id]["_".join(species_homology.split("_")[0:2]).capitalize()] != 0:
            if species_homology == species:
                continue
            common_ancestor1 = repr(tree.common_ancestor("_".join(species_homology.split("_")[0:2]).capitalize(),
                                                         "_".join(species.split("_")[0:2]).capitalize()))
            if 'Clade()' in common_ancestor1:
                if trace_main_branch[common_ancestor1] != 0:
                    for step in all_trace_dic["_".join(species_homology.split("_")[0:2]).capitalize()]:
                        gene_homology.append(step)

            elif 'name' in common_ancestor1:
                if trace_main_branch[common_ancestor1.split("name='")[1].split("')")[0]] != 0:
                    for step in all_trace_dic["_".join(species_homology.split("_")[0:2]).capitalize()]:
                        gene_homology.append(step)
            else:
                if trace_main_branch['confidence=' + common_ancestor1.split("confidence=")[1].split(")")[0]] != 0:
                    for step in all_trace_dic["_".join(species_homology.split("_")[0:2]).capitalize()]:
                        gene_homology.append(step)

    if ori_branch == 'Clade()':
        ori_branch = 'confidence=None'

    if ori_branch == 'confidence=None':
        gene_age = ">" + str(maximum_age)
        gene_interval = ">" + str(maximum_age)
        gene_branch = 0
    else:
        gene_branch = len(trace) - trace.index(ori_branch) - 1

        age_lower = 0
        for branch in trace[0:trace.index(ori_branch)]:
            age_lower += branch_length_dict[branch]
        age_lower = int(age_lower)
        age_upper = 0
        for branch in trace[0:trace.index(ori_branch) + 1]:
            age_upper += branch_length_dict[branch]
        age_upper = int(age_upper)
        gene_age = (age_lower + age_upper) / 2

        gene_interval = str(age_lower) + "-" + str(age_upper)
    print(gene_age, gene_interval, gene_id)
    return gene_age, gene_interval, gene_branch, ori_branch, gene_homology


def run_count(homologyB, nwk_dir, param, species, division):
    # Call as java ca.umontreal.iro.evolution.genecontent.AsymmetricWagner[-gain x] [-max_paralogs n][-min_lineages]
    # phylogeny table]
    # homology B:species:gene ID:homology species

    family2node2presence = {}
    species_list = [" ".join(i.split("_")[0:2]).capitalize() for i in
                    homologyB[species][list(homologyB[species].keys())[0]]]
    column = "Family" + "\\t" + "\\t".join(
        [" ".join(i.split("_")[0:2]).capitalize() for i in species_list])
    times = ""
    for gene_id in homologyB[species]:
        homology_list = []
        times_line = [gene_id]
        for homology_len in species_list:
            homology_len = "_".join(homology_len.split(" "))
            if homology_len not in homologyB[species][gene_id]:
                times_line.append("0")
            else:
                times_line.append(str(len(set(homologyB[species][gene_id][homology_len]))))

        # homologyB[species][gene_id] = homology_list
        times_line = "\\t".join(times_line)
        times_line = times_line + "\\n"
        times = times + times_line
    times.strip("\\n")
    times.strip("\\t")
    result = column.strip("\\t") + "\\n" + times.strip("\\t")
    temp_dir = work_dir + "temp_cookie/111" + species + "." + gene_id + "." + division + ".tsv"
    file_content = '''%s''' % result
    outfile = open(temp_dir, "w")
    outfile.write(file_content)
    outfile.close()

    # param = " ".join(sys.argv[1:])
    if param == "":
        print("Gain penalty = -gain 5.0")
    else:
        print(species, "Gain penalty : " + param)
    cmd = 'java -Xmx2048M -cp /home/lchen/genorigin/Count/Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner %s %s %s' % (
        param, nwk_dir, temp_dir)
    output = subprocess.getoutput(cmd)
    nodes = []

    for line in output.split('\\n'):
        if "FAMILY" in line:
            t = line.split('\\t')
            if t[1] == 'name':
                nodes = [x.replace(' ', '_') for x in t[2:-4]]
                nodes[-1] = 'Clade()'
            else:
                counts = [int(x) for x in t[2:-4]]
                family2node2presence[t[1]] = dict(zip(nodes, counts))
    os.remove(temp_dir)

    return family2node2presence


def common_ancestor_branch(other, tree):
    try:
        branch = tree.get_path("_".join(species.split("_")[0:2]).capitalize()).index(
            tree.common_ancestor("_".join(species.split("_")[0:2]).capitalize(),
                                 "_".join(other.strip().split("_")[0:2]).capitalize())) + 1
    except:
        branch = 0
    return branch


# through count analyse gene age
def main(homologyB, nwk_dir, family2node2presence, species, species_list, division):
    print(species, 'start main')
    tree = Phylo.read(nwk_dir, "newick")
    species_tree_path = tree.get_path("_".join(species.split("_")[0:2]).capitalize())

    # branch_species_dict
    all_trace_dic = {}
    branch_species_dict = defaultdict(lambda: [])
    # a = species_trace,
    for species_use in species_list:
        species_use_name = "_".join(species_use.split("_")[0:2]).capitalize()
        branch_species_dict[common_ancestor_branch(species_use_name, tree)].append(species_use_name)
        a = [j.split(")")[0] for j in [i.split(", ")[1] for i in
                                       repr(tree.get_path(species_use_name)).split(
                                           "Clade(")[1:]]]
        a[-1] = a[-1].split("'")[1]
        all_trace_dic[species_use_name] = a

    '''try:
        max_mid = repr(species_tree_path).split("Clade")[1:]
        max_mid.reverse()
        max_branch = repr(tree.common_ancestor(["_".join(i.split("_")[0:2]).capitalize() for i in os.listdir(
            work_dir + "homology_json/" + division)])).split("confidence=")[1].split(")")[0]
        max_list = []
        for i in max_mid:
            max_list.append(i.split("branch_length=")[1].split(",")[0])
            if "confidence=" + max_branch in i:
                break
        maximum_age = sum([float(i) for i in max_list])
    except:
        maximum_age = 4290'''
    maximum_age = 4290
    branch_length = [float(i.split(",")[0]) for i in repr(species_tree_path).split('branch_length=')[1:]]

    trace = all_trace_dic["_".join(species.split("_")[0:2]).capitalize()]

    branch_length_dict = dict(zip(trace, branch_length))

    trace.reverse()
    trace.append("Clade()")

    # analyse_group_max_age
    for group in analyse_group[division]:
        if "_".join(species.split("_")[0:2]).capitalize() in group:
            analyse_group_root = repr(tree.common_ancestor(group))
    if 'confidence=' in analyse_group_root:
        analyse_group_root = analyse_group_root.split("confidence=")[1].split(")")[0]
    elif "name" in analyse_group_root:
        analyse_group_root = analyse_group_root.split("name='")[1].split("')")[0]
    elif "Clade()" in analyse_group_root:
        analyse_group_root = "Clade()"

    species_path = repr(tree.get_path("_".join(species.split("_")[0:2]).capitalize())).split("Clade(")[1:]
    aa = []
    for i in species_path:
        if 'confidence=' in i:
            aa.append(i.split("confidence=")[1].split(")")[0])
        elif "name" in i:
            aa.append(i.split("name='")[1].split("')")[0])
    species_path = aa
    for node in species_path:
        if analyse_group_root == node:
            analyse_group_root_node_index = species_path.index(node)
            analyse_group_max_age_list = []
            for i in species_path[analyse_group_root_node_index + 1:]:
                if i.isdigit():
                    analyse_group_max_age_list.append(branch_length_dict["confidence=" + str(i)])
                else:
                    analyse_group_max_age_list.append(branch_length_dict[i])
            count = 0
            for i in analyse_group_max_age_list:
                count += float(i)
            analyse_group_max_age = count

    print('species2pan_homology')
    # species2pan_homology
    if species not in pan_species_list:
        species2pan_homology = read_json(work_dir + 'species_to_pan_species/' + species)

        pan_list = []
        for i in species2pan_homology:
            for j in species2pan_homology[i]:
                for k in species2pan_homology[i][j]:
                    for l in species2pan_homology[i][j][k]:
                        pan_list.append(species2species[species2pan_homology[i][j][k][l][0]])
        pan_list = list(set(pan_list))
        pan_homology = {}
        for pan_species in pan_list:
            pan_homology[pan_species] = read_json(work_dir + 'pan_extension/' + pan_species)[pan_species]
    else:
        species2pan_homology = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: []))))
        for gene_id in homologyB[species]:
            species2pan_homology["_".join(species.split("_")[0:2]).capitalize()][gene_id]['100'][gene_id] = [
                "_".join(species.split("_")[0:2]).capitalize()]

        pan_homology = {species: read_json(work_dir + 'pan_extension/' + species)[species]}

    # function
    function = defaultdict(lambda: defaultdict(str))

    ori_branch = "_".join(species.split("_")[0:2]).capitalize()

    for gene_id in homologyB[species]:
        mid = dict(zip(total_species_list, [[] for i in total_species_list]))
        mid.update(homologyB[species][gene_id])
        homologyB[species][gene_id] = mid
        effect = []
        for node in family2node2presence[gene_id]:
            if family2node2presence[gene_id][node] > 0:
                if node.isdigit():
                    effect.append('confidence=' + str(node))
                else:
                    effect.append(node)
        trace_main_branch = {}
        for step in trace:
            if step in effect:
                ori_branch = step
                trace_main_branch[step] = 1
            else:
                break

        gene_homology = []
        for step in trace:
            gene_homology.append(step)
            if step not in trace_main_branch:
                trace_main_branch[step] = 0

        for species_homology in species_list:
            if family2node2presence[gene_id]["_".join(species_homology.split("_")[0:2]).capitalize()] != 0:
                if species_homology == species:
                    continue
                common_ancestor1 = repr(tree.common_ancestor("_".join(species_homology.split("_")[0:2]).capitalize(),
                                                             "_".join(species.split("_")[0:2]).capitalize()))
                if 'Clade()' in common_ancestor1:
                    if trace_main_branch[common_ancestor1] != 0:
                        for step in all_trace_dic["_".join(species_homology.split("_")[0:2]).capitalize()]:
                            gene_homology.append(step)

                elif 'name' in common_ancestor1:
                    if trace_main_branch[common_ancestor1.split("name='")[1].split("')")[0]] != 0:
                        for step in all_trace_dic["_".join(species_homology.split("_")[0:2]).capitalize()]:
                            gene_homology.append(step)
                else:
                    if trace_main_branch['confidence=' + common_ancestor1.split("confidence=")[1].split(")")[0]] != 0:
                        for step in all_trace_dic["_".join(species_homology.split("_")[0:2]).capitalize()]:
                            gene_homology.append(step)

        if ori_branch == 'Clade()':
            ori_branch = 'confidence=None'

        if ori_branch == 'confidence=None':
            gene_age = ">" + str(maximum_age)
            gene_interval = ">" + str(maximum_age)
            gene_branch = 0
        else:
            gene_branch = len(trace) - trace.index(ori_branch) - 1

            age_lower = 0
            for branch in trace[0:trace.index(ori_branch)]:
                age_lower += branch_length_dict[branch]
            age_lower = int(age_lower)
            age_upper = 0
            for branch in trace[0:trace.index(ori_branch) + 1]:
                age_upper += branch_length_dict[branch]
            age_upper = int(age_upper)
            gene_age = (age_lower + age_upper) / 2

            gene_interval = str(age_lower) + "-" + str(age_upper)

        mid = gene_age
        if type(gene_age) == float:
            pass
        else:
            gene_age = gene_age.split(">")[1]

        if float(gene_age) >= float(analyse_group_max_age):
            gene_age = mid
            (gene_age, gene_interval, gene_branch, ori_branch, gene_homology, homology_changed
             ) = total_tree(gene_id, species, division, nwk_dir, trace, branch_length_dict, maximum_age,
                            all_trace_dic, gene_age, gene_interval, gene_branch, ori_branch, gene_homology,
                            homologyB[species][gene_id], species2pan_homology, pan_homology)
            homologyB[species][gene_id] = homology_changed

            function[gene_id]["ensembl_gene_id"] = gene_id
            function[gene_id]["scientific_name"] = species
            if gene_branch == 0:
                gene_age = ">4290"
            function[gene_id]["gene_age"] = gene_age
            function[gene_id]["gene_branch"] = gene_branch
            function[gene_id]["gene_homolog"] = list(set(gene_homology))
            function[gene_id]["ori_branch"] = ori_branch
            function[gene_id]["trace"] = trace
            function[gene_id]["gene_Interval"] = gene_interval

        else:
            function[gene_id]["ensembl_gene_id"] = gene_id
            function[gene_id]["scientific_name"] = species
            function[gene_id]["gene_age"] = gene_age
            function[gene_id]["gene_branch"] = gene_branch
            function[gene_id]["gene_homolog"] = list(set(gene_homology))
            function[gene_id]["ori_branch"] = ori_branch
            function[gene_id]["trace"] = trace
            function[gene_id]["gene_Interval"] = gene_interval

    function_dir = work_dir + "mid_merge/" + species + "_" + gene_id + "_homology"
    file_content = '''%s''' % json.dumps(homologyB)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    function_dir = work_dir + "mid_merge/" + species + "_" + gene_id + "_function"
    file_content = '''%s''' % json.dumps(function)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()
    print(species, "end")


def run(hh):
    division, species, gene_list = hh[0], hh[1], hh[2]

    species_list = list(set(["_".join(i.split("_")[0:2]).capitalize() for j in
                             os.listdir(work_dir + "homology_json/") for i
                             in os.listdir(work_dir + "homology_json/" + j)]))
    global homologyB
    homologyB_dir = work_dir + "homology_json/" + division + "/" + species

    # homologyB = read_json(homologyB_dir)

    homologyB = """ + json.dumps({species: dict(zip(list(homologyB[species].keys())[i * 500:(i + 1) * 500],
                                                               list(homologyB[species].values())[
                                                               i * 500:(i + 1) * 500]))})

    python_file2 = """
    family2node2presence = run_count(homologyB, nwk_dir, param, species, division)

    main(homologyB, nwk_dir, family2node2presence, species, species_list, division)

species2species = {}
for division in os.listdir(work_dir + "homology_json/"):
    for species in os.listdir(work_dir + "homology_json/" + division):
        species2species['_'.join(species.split('_')[0:2]).capitalize()] = species

param = "-gain 1.4"

nwk_dir = work_dir + "total.nwk"

species = '%s'

division = '%s'

total_species_list = list(set(["_".join(i.split("_")[0:2]).capitalize() for j in
                               os.listdir(work_dir + "homology_json/") for i
                               in os.listdir(work_dir + "homology_json/" + j)]))
pan_species_list = os.listdir(work_dir + "homology_json/pan_homology")

run([division, species, []])
    """ % (species, division)
    function_dir = work_dir + 'python_all_homology/'+list(homologyB[species].keys())[i]+'.py'
    file_content = python_file1+python_file2
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()
