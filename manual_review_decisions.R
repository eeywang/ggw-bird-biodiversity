# =============================================================================
# manual_review_decisions.R — Hand-curated species decisions for 3 focal countries
#
# These decisions were made by manually checking each species' ecology,
# range, and taxonomy. They are applied in 04_apply_review.R.
#
# Decision codes:
#   "remove"         — confirmed non-African / extinct / genus-level / pelagic
#   "keep_resident"  — confirmed African resident (sedentary)
#   "keep_partial"   — confirmed partial migrant or Palearctic winterer
# =============================================================================

# ============================================================================
# NIGERIA
# ============================================================================

# Nigeria went through iterative V2→V3→V4 reviews. The decisions below
# consolidate the final V4 outcome.

# Unknown species that were resolved in V3/V4 review rounds:
nigeria_unknown_decisions <- data.frame(
  species = c(
    # V4: unknown species resolved
    "Phylloscopus sibillatrix",   # misspelling → remove
    "Chroicocephalus saundersi",  # non-African → remove
    "Vanellus macropterus",       # extinct → remove
    "Tadorna cristata",           # extinct → remove
    "Haematopus meadewaldoi",     # extinct → remove
    "Rhodonessa caryophyllacea",  # extinct → remove
    "Prosobonia cancellata",      # extinct → remove
    "Prosobonia leucoptera",      # extinct → remove
    "Leptoptilos crumeniferus",   # African resident
    "Rhinoptilus africanus",      # African resident
    "Riparia cincta",             # African resident
    "Pseudalethe poliocephala",   # African resident
    "Pseudoalcippe abyssinica",   # African resident
    "Salpornis spilonotus",       # African resident
    "Anthus longicaudatus",       # African resident
    "Petrochelidon fuliginosa"    # African resident
  ),
  decision = c(
    "remove", "remove", "remove", "remove",
    "remove", "remove", "remove", "remove",
    "keep_resident", "keep_resident", "keep_resident",
    "keep_resident", "keep_resident", "keep_resident",
    "keep_resident", "keep_resident"
  ),
  stringsAsFactors = FALSE
)

# Full migrants to keep (have African resident populations or are
# Palearctic winterers already captured in shared_lists$palearctic_winterers)
nigeria_migrants_to_keep <- c(
  # V3 fix: species erroneously removed as full migrants
  "Crex egregia",             # African Crake — African resident
  "Porphyrio alleni",         # Allen's Gallinule — African resident
  "Milvus aegyptius",         # Yellow-billed Kite — African resident/partial
  "Aquila rapax",             # Tawny Eagle — African resident population
  "Ardea purpurea",           # Purple Heron — African resident population
  "Ardeola ralloides",        # Squacco Heron — African resident population
  "Nycticorax nycticorax",    # Night Heron — African resident population
  "Ixobrychus minutus",       # Little Bittern — African resident population
  # V3 fix: Palearctic winterers missed in initial coding
  "Cercotrichas galactotes",  # = Erythropygia galactotes
  "Curruca communis",         # Common Whitethroat (= Sylvia communis)
  "Curruca curruca",          # Lesser Whitethroat (= Sylvia curruca)
  "Phylloscopus sibilatrix",  # Wood Warbler (correct spelling)
  "Crex crex"                 # Corncrake — Palearctic winterer
)

# Erroneous records (non-African species in GBIF for Nigeria)
nigeria_error_records <- c()


# ============================================================================
# SENEGAL
# ============================================================================

senegal_unknown_decisions <- data.frame(
  species = c(
    NA,                          # invalid record
    "Anarhynchus pecuarius",     # Kittlitz's Plover
    "Erythropygia galactotes",   # Rufous-tailed Scrub Robin
    "Thinornis dubius",          # = Charadrius pecuarius
    "Tauraco violaceus",         # Violet Turaco
    "Anarhynchus alexandrinus",  # Kentish Plover
    "Cecropis rufula",           # Red-rumped Swallow
    "Turdoides fulva",           # Fulvous Babbler
    "Puffinus griseus",          # Sooty Shearwater — pelagic
    "Pentholaea albifrons",      # White-fronted Black Chat
    "Puffinus gravis",           # Great Shearwater — pelagic
    "Aerospiza tachiro",         # African Goshawk
    "Amirafra rufocinnamomea",   # Rufous-naped Lark
    "Oceanodroma castro",        # Band-rumped Storm Petrel — pelagic
    "Egretta dimorpha",          # Dimorphic Egret
    "Oceanodroma leucorhoa",     # Leach's Storm Petrel — pelagic
    "Sterna"                     # genus-level record
  ),
  decision = c(
    "remove",           # NA
    "keep_resident",    # Anarhynchus pecuarius
    "keep_partial",     # Erythropygia galactotes
    "keep_resident",    # Thinornis dubius
    "keep_resident",    # Tauraco violaceus
    "keep_partial",     # Anarhynchus alexandrinus
    "keep_resident",    # Cecropis rufula
    "keep_resident",    # Turdoides fulva
    "remove",           # Puffinus griseus (pelagic)
    "keep_resident",    # Pentholaea albifrons
    "remove",           # Puffinus gravis (pelagic)
    "keep_resident",    # Aerospiza tachiro
    "keep_resident",    # Amirafra rufocinnamomea
    "remove",           # Oceanodroma castro (pelagic)
    "keep_resident",    # Egretta dimorpha
    "remove",           # Oceanodroma leucorhoa (pelagic)
    "remove"            # Sterna (genus-level)
  ),
  stringsAsFactors = FALSE
)

senegal_migrants_to_keep <- c(
  # African resident populations
  "Milvus migrans",             # Black Kite — African subspecies
  "Apus affinis",               # Little Swift
  "Cisticola juncidis",         # Zitting Cisticola
  "Tachyspiza badia",           # Shikra
  "Plectropterus gambensis",    # Spur-winged Goose
  "Botaurus minutus",           # Little Bittern (= Ixobrychus minutus)
  "Accipiter badius",           # Shikra
  "Upupa africana",             # African Hoopoe
  "Acrocephalus baeticatus",    # African Reed Warbler
  "Apus caffer",                # White-rumped Swift
  # Intra-African migrants / partial
  "Pandion haliaetus",          # Osprey
  "Falco tinnunculus",          # Common Kestrel
  "Lanius excubitor",           # Great Grey Shrike
  "Charadrius dubius",          # Little Ringed Plover
  "Cinnyricinclus leucogaster", # Violet-backed Starling
  "Clamator jacobinus",         # Jacobin Cuckoo
  "Clamator levaillantii",      # Levaillant's Cuckoo
  "Clamator glandarius",        # Great Spotted Cuckoo
  "Glareola pratincola",        # Collared Pratincole
  "Iduna pallida"               # Eastern Olivaceous Warbler
)

senegal_error_records <- c(
  "Haliaeetus leucocephalus",   # Bald Eagle — North American!
  "Fregata magnificens",        # Magnificent Frigatebird — Americas!
  "Spatula discors"             # Blue-winged Teal — North American
)


# ============================================================================
# ETHIOPIA
# ============================================================================

ethiopia_unknown_decisions <- data.frame(
  species = c(
    NA,                          # invalid record
    "Crinifer leucogaster",      # White-bellied Go-away-bird
    "Buphagus erythroryncha",    # Red-billed Oxpecker
    "Turdus simensis",           # Ethiopian Thrush (endemic)
    "Dioptrornis chocolatinus",  # Abyssinian Slaty Flycatcher
    "Turdoides rubiginosa",      # Rufous Chatterer
    "Erythropygia leucophrys",   # White-browed Scrub Robin
    "Menelikornis leucotis",     # Northern Pied Hornbill
    "Pogoniulus uropygialis",    # White-rumped Tinkerbird
    "Parus thruppi",             # Somali Tit
    "Crinifer personatus",       # Bare-faced Go-away-bird
    "Tachyspiza minulla",        # African Little Sparrowhawk
    "Parus leuconotus",          # White-backed Tit (endemic)
    "Aerospiza tachiro",         # African Goshawk
    "Erythropygia galactotes",   # Rufous-tailed Scrub Robin
    "Menelikornis ruspolii",     # Ruspoli's Hornbill (endemic)
    "Anarhynchus pecuarius",     # Kittlitz's Plover
    "Thinornis dubius",          # Forbes's Plover
    "Calendulauda gilletti",     # Gillett's Lark
    "Grus carunculatus",         # Wattled Crane
    "Turdoides aylmeri",         # Scaly Chatterer
    "Calandrella somalica",      # Somali Short-toed Lark
    "Pseudalaemon fremantlii",   # Short-tailed Lark
    "Poeoptera sharpii",         # Sharpe's Starling
    "Calandrella erlangeri",     # Erlanger's Lark (endemic)
    "Amirafra rufocinnamomea",   # Rufous-naped Lark
    "Anarhynchus alexandrinus",  # Kentish Plover
    "Pentholaea albifrons",      # White-fronted Black Chat
    "Corypha hypermetra",        # Red-winged Lark
    "Anarhynchus asiaticus",     # Caspian Plover
    "Anarhynchus atrifrons",     # Black-fronted Dotterel
    "Empidornis semipartitus",   # Silverbird
    "Anarhynchus leschenaultii", # Greater Sand Plover
    "Cecropis rufula",           # Red-rumped Swallow
    "Phoeniculus granti",        # Grant's Wood Hoopoe
    "Corypha kidepoensis",      # Kidepo Lark
    "Tachyspiza brevipes"        # Ovambo Sparrowhawk
  ),
  decision = c(
    "remove",           # NA
    "keep_resident",    # Crinifer leucogaster
    "keep_resident",    # Buphagus erythroryncha
    "keep_resident",    # Turdus simensis
    "keep_resident",    # Dioptrornis chocolatinus
    "keep_resident",    # Turdoides rubiginosa
    "keep_resident",    # Erythropygia leucophrys
    "keep_resident",    # Menelikornis leucotis
    "keep_resident",    # Pogoniulus uropygialis
    "keep_resident",    # Parus thruppi
    "keep_resident",    # Crinifer personatus
    "keep_resident",    # Tachyspiza minulla
    "keep_resident",    # Parus leuconotus
    "keep_resident",    # Aerospiza tachiro
    "keep_partial",     # Erythropygia galactotes
    "keep_resident",    # Menelikornis ruspolii
    "keep_resident",    # Anarhynchus pecuarius
    "keep_resident",    # Thinornis dubius
    "keep_resident",    # Calendulauda gilletti
    "keep_resident",    # Grus carunculatus
    "keep_resident",    # Turdoides aylmeri
    "keep_resident",    # Calandrella somalica
    "keep_resident",    # Pseudalaemon fremantlii
    "keep_resident",    # Poeoptera sharpii
    "keep_resident",    # Calandrella erlangeri
    "keep_resident",    # Amirafra rufocinnamomea
    "keep_partial",     # Anarhynchus alexandrinus
    "keep_resident",    # Pentholaea albifrons
    "keep_resident",    # Corypha hypermetra
    "keep_partial",     # Anarhynchus asiaticus
    "keep_resident",    # Anarhynchus atrifrons
    "keep_resident",    # Empidornis semipartitus
    "keep_partial",     # Anarhynchus leschenaultii
    "keep_resident",    # Cecropis rufula
    "keep_resident",    # Phoeniculus granti
    "keep_resident",    # Corypha kidepoensis
    "keep_resident"     # Tachyspiza brevipes
  ),
  stringsAsFactors = FALSE
)

ethiopia_migrants_to_keep <- c(
  # African resident populations
  "Milvus migrans",             # Black Kite
  "Falco tinnunculus",          # Common Kestrel
  "Neophron percnopterus",      # Egyptian Vulture
  "Anthus cinnamomeus",         # African Pipit — AVONET error!
  "Lamprotornis shelleyi",      # Shelley's Starling — AVONET error!
  "Tachyspiza badia",           # Shikra
  "Plectropterus gambensis",    # Spur-winged Goose
  "Tachymarptis melba",         # Alpine Swift
  "Apus affinis",               # Little Swift
  "Apus caffer",                # White-rumped Swift
  "Cisticola juncidis",         # Zitting Cisticola
  "Botaurus minutus",           # Little Bittern
  "Accipiter badius",           # Shikra
  "Acrocephalus baeticatus",    # African Reed Warbler
  "Upupa africana",             # African Hoopoe
  "Psalidoprocne albiceps",     # White-headed Saw-wing
  "Petrochelidon spilodera",    # South African Cliff Swallow
  # Intra-African migrants / winterers
  "Iduna pallida",              # Eastern Olivaceous Warbler
  "Lanius excubitor",           # Great Grey Shrike
  "Cinnyricinclus leucogaster", # Violet-backed Starling
  "Clamator jacobinus",         # Jacobin Cuckoo
  "Clamator levaillantii",      # Levaillant's Cuckoo
  "Clamator glandarius",        # Great Spotted Cuckoo
  "Glareola pratincola",        # Collared Pratincole
  "Cuculus clamosus",           # Black Cuckoo
  "Coturnix delegorguei"        # Harlequin Quail
)

ethiopia_error_records <- c()
# No obvious non-African error records found for Ethiopia

# ---- Bundle all decisions ----
country_decisions <- list(
  Nigeria = list(
    unknown_decisions  = nigeria_unknown_decisions,
    migrants_to_keep   = nigeria_migrants_to_keep,
    error_records      = nigeria_error_records
  ),
  Senegal = list(
    unknown_decisions  = senegal_unknown_decisions,
    migrants_to_keep   = senegal_migrants_to_keep,
    error_records      = senegal_error_records
  ),
  Ethiopia = list(
    unknown_decisions  = ethiopia_unknown_decisions,
    migrants_to_keep   = ethiopia_migrants_to_keep,
    error_records      = ethiopia_error_records
  )
)

cat("=== manual_review_decisions.R loaded ===\n")
cat("  Decisions defined for:", paste(names(country_decisions), collapse = ", "), "\n")
