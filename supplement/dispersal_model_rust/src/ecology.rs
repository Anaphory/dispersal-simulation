use crate::KCal;

use std::collections::HashMap;

pub fn load_precipitation_tif() -> Option<(Vec<u16>, u32)> {
    // TODO: How do I specify data paths?
    let f = std::fs::File::open("wc2.1_5m_bio_12-16bit.tif").ok()?;
    let maybe_image = tiff::decoder::Decoder::new(f);
    let mut image = maybe_image.ok()?;
    let (width, height) = image.dimensions().ok()?;
    println!("# {}x{}", width, height);
    let outcome = image.read_image().ok()?;
    println!("# Image read");
    let vec = match outcome {
        tiff::decoder::DecodingResult::U8(v) => v.iter().map(|g| u16::from(*g)).collect(),
        tiff::decoder::DecodingResult::U16(w) => w,
        tiff::decoder::DecodingResult::U32(v) => v.iter().map(|g| *g as u16).collect(),
        tiff::decoder::DecodingResult::U64(v) => v.iter().map(|g| *g as u16).collect(),
    };
    Some((vec, width))
}

pub fn load_density_tif() -> Option<(Vec<u16>, u32)> {
    let f = std::fs::File::open("../tallavaara/lpdensity.tif").ok()?;
    let maybe_image = tiff::decoder::Decoder::new(f);
    let mut image = maybe_image.ok()?;
    let (width, height) = image.dimensions().ok()?;
    println!("# {}x{}", width, height);
    let outcome = image.read_image().ok()?;
    println!("# Image read");
    let vec = match outcome {
        tiff::decoder::DecodingResult::U8(v) => v.iter().map(|g| u16::from(*g)).collect(),
        tiff::decoder::DecodingResult::U16(w) => w,
        tiff::decoder::DecodingResult::U32(v) => v.iter().map(|g| *g as u16).collect(),
        tiff::decoder::DecodingResult::U64(v) => v.iter().map(|g| *g as u16).collect(),
    };
    Some((vec, width))
}

pub const ATTESTED_ECOREGIONS: usize = 303;

#[allow(clippy::excessive_precision, clippy::unreadable_literal)]
pub fn patch_from_ecoregions(ecoregion: i64, area_in_100_km2: f32) -> KCal {
    let logpopdensity: HashMap<_, KCal> = vec![
        (2, 4.0974917821114),       // Admiralty Islands lowland rain forests
        (3, 2.13331596527082),      // Aegean and Western Turkey sclerophyllous and mixed forests
        (4, 1.30248547521061),      // Afghan Mountains semi-desert
        (5, 0.977186366941075),     // Ahklun and Kilbuck Upland Tundra
        (6, 0.777874727210317),     // Al-Hajar foothill xeric woodlands and shrublands
        (7, 0.759580174873787),     // Al-Hajar montane woodlands and shrublands
        (8, 1.6117677865496),       // Alai-Western Tian Shan steppe
        (9, 0.463267443958331),     // Alashan Plateau semi-desert
        (10, 0.633686524205738),    // Alaska-St. Elias Range tundra
        (11, 1.71320577037012),     // Alaska Peninsula montane taiga
        (12, 2.95334229676252),     // Albany thickets
        (13, 1.86421750985489),     // Alberta-British Columbia foothills forests
        (14, 1.47547851157093),     // Albertine Rift montane forests
        (15, 1.08302391180798),     // Aldabra Island xeric scrub
        (16, 1.84050212615075),     // Aleutian Islands tundra
        (17, 3.14616602870839),     // Allegheny Highlands forests
        (18, 2.46642781775555),     // Alps conifer and mixed forests
        (19, 1.71548985598458),     // Altai alpine meadow and tundra
        (20, 1.4546657063152),      // Altai montane forest and forest steppe
        (21, 1.94369241997949),     // Altai steppe and semi-desert
        (22, 2.96502991730754),     // Alto Paraná Atlantic forests
        (23, 2.95309228707746),     // Amazon-Orinoco-Southern Caribbean mangroves
        (25, 1.91768844789713),     // Amur meadow steppe
        (26, 2.1610548492277),      // Anatolian conifer and deciduous mixed forests
        (27, 2.74517811949172),     // Andaman Islands rain forests
        (28, 1.89044088368986),     // Angolan montane forest-grassland
        (29, 1.94002555730349),     // Angolan mopane woodlands
        (30, 1.59672727070173),     // Angolan scarp savanna and woodlands
        (31, 1.55318109302441),     // Angolan wet miombo woodlands
        (33, 3.46165336032365),     // Appalachian-Blue Ridge forests
        (34, 3.52347164518471),     // Appalachian mixed mesophytic forests
        (35, 3.39146333337869),     // Appalachian Piedmont forests
        (36, 2.09104132905453),     // Appenine deciduous montane forests
        (37, 2.70759930310593),     // Apure-Villavicencio dry forests
        (38, 0.422137324741384),    // Arabian sand desert
        (39, 2.91306273146247),     // Araucaria moist forests
        (40, 2.95907106025888),     // Araya and Paria xeric scrub
        (41, -0.416250895572695),   // Arctic coastal tundra
        (42, -0.482022188667833),   // Russian Arctic desert
        (43, -0.324673473197806),   // Arctic foothills tundra
        (44, 2.31608968883798),     // Arizona Mountains forests
        (45, 3.01581550533565),     // Arnhem Land tropical savanna
        (47, 0.23848318791909),     // Atacama desert
        (48, 2.7889888881731),      // Atlantic Coast restingas
        (49, 0.650535863983308),    // Saharan Atlantic coastal desert
        (50, 3.87323238446442),     // Atlantic coastal pine barrens
        (51, 0.0750922889069491),   // Congolian coastal forests
        (52, 2.41772719851727),     // European Atlantic mixed forests
        (53, 3.09091638124688),     // Australian Alps montane grasslands
        (54, 2.06903492885911),     // Azerbaijan shrub desert and steppe
        (55, 2.01095505938591),     // Azores temperate mixed forests
        (56, 1.36647861006464),     // Badghyz and Karabil semi-desert
        (57, 2.53265332466742),     // Bahamian-Antillean mangroves
        (58, 2.96831159345687),     // Bahia coastal forests
        (59, 2.78429291591059),     // Bahia interior forests
        (60, 0.811593287213633),    // Baja California desert
        (61, 2.73736464555287),     // Bajío dry forests
        (62, 2.00917525364037),     // Balkan mixed forests
        (63, 2.63089248721969),     // Balsas dry forests
        (64, 2.26086191368584),     // Baltic mixed forests
        (65, 1.25976866468217),     // Baluchistan xeric woodlands
        (66, 3.36167042980405),     // Banda Sea Islands moist deciduous forests
        (67, 2.71167022610384),     // Belizian pine savannas
        (68, 3.48534866747345),     // Beni savanna
        (69, 0.196243263610247),    // Russian Bering tundra
        (70, 1.01465370074617),     // Beringia lowland tundra
        (71, 0.326457896489572),    // Beringia upland tundra
        (73, 4.2150743753271),      // Biak-Numfoor rain forests
        (74, 2.52163101562382),     // Blue Mountains forests
        (75, 1.45675495681477),     // Bohai Sea saline meadow
        (76, 2.5959591173798),      // Bolivian montane dry forests
        (77, 3.69126178064036),     // Bolivian Yungas
        (78, 2.51607108095098),     // Borneo lowland rain forests
        (79, 2.80020374724826),     // Borneo montane rain forests
        (80, 2.36078731622686),     // Borneo peat swamp forests
        (81, 0.578821318460697),    // Brahmaputra Valley semi-evergreen forests
        (82, 2.5661359466275),      // Brigalow tropical savanna
        (83, 2.69419026108462),     // British Columbia coastal conifer forests
        (84, -0.287086734179644),   // Brooks-British Range tundra
        (85, 2.90179776960045),     // Buru rain forests
        (86, 2.67612682978226),     // Caatinga
        (87, 2.55248957277183),     // Caatinga Enclaves moist forests
        (88, 2.3386809590726),      // Caledon conifer forests
        (89, 1.9583882329741),      // California Central Valley grasslands
        (90, 1.92474059543653),     // California coastal sage and chaparral
        (91, 2.87221825678334),     // California interior chaparral and woodlands
        (92, 2.41399338782549),     // California montane chaparral and woodlands
        (93, -0.0381918068385647),  // Cameroon Highlands forests
        (94, 2.56139208621003),     // Campos Rupestres montane savanna
        (95, 2.1333310983066),      // Canadian Aspen forests and parklands
        (96, -0.460116747381652),   // Canadian Low Arctic tundra
        (97, 0.472555074132789),    // Canary Islands dry woodlands and forests
        (98, 2.80777521615694),     // Cantabrian mixed forests
        (99, 2.5225830487383),      // Canterbury-Otago tussock grasslands
        (100, 0.566521564058258),   // Cape Verde Islands dry forests
        (101, 3.34053202865537),    // Cape York Peninsula tropical savanna
        (102, 2.51653114137904),    // Caqueta moist forests
        (103, 1.79810499746683),    // Cardamom Mountains rain forests
        (104, 2.16781252622051),    // Caribbean shrublands
        (105, 1.03758617843162),    // Carnarvon xeric shrublands
        (106, 3.98194780349969),    // Carolines tropical moist forests
        (107, 2.21976928082319),    // Carpathian montane forests
        (108, 2.35232661638639),    // Carpentaria tropical savanna
        (109, 2.2624539673653),     // Caspian Hyrcanian mixed forests
        (110, 1.00422146384796),    // Caspian lowland desert
        (111, 2.83000858307821),    // Catatumbo moist forests
        (112, 2.78885833196344),    // Cauca Valley dry forests
        (113, 2.81075698613224),    // Cauca Valley montane forests
        (114, 2.5358515540056),     // Caucasus mixed forests
        (116, 2.52791891576037),    // Celtic broadleaf forests
        (117, 1.09935829945159),    // Central Afghan Mountains xeric woodlands
        (118, -0.75501333355935),   // Central African mangroves
        (119, 2.83417087977336),    // Central American Atlantic moist forests
        (120, 2.64673274131715),    // Central American dry forests
        (121, 2.91083070199516),    // Central American montane forests
        (122, 2.77189408722518),    // Central American pine-oak forests
        (123, 1.75972247999263),    // Central Anatolian steppe
        (124, 1.89663044541613),    // Central Anatolian steppe and woodlands
        (125, 3.53104022389716),    // Central-Southern Cascades Forests
        (126, 2.71150295723863),    // Central-Southern US mixed grasslands
        (127, 0.995947415823361),   // Central Andean dry puna
        (128, 1.68040671224525),    // Central Andean puna
        (129, 2.61204580625724),    // Central Andean wet puna
        (130, 0.887190712687629),   // Central Asian northern desert
        (131, 0.996414459754437),   // Central Asian riparian woodlands
        (132, 0.837547681341407),   // Central Asian southern desert
        (133, 1.8589780409135),     // Central British Columbia Mountain forests
        (134, 2.58130324703487),    // Central bushveld
        (135, 1.55892034992629),    // Central Canadian Shield forests
        (136, 1.6263965527289),     // Central China Loess Plateau mixed forests
        (137, 0.230600513408591),   // Central Congolian lowland forests
        (138, 0.886799521223265),   // Central Deccan Plateau dry deciduous forests
        (139, 2.22833654521916),    // Central European mixed forests
        (140, 3.550187927864),      // Central US forest-grasslands transition
        (141, 1.61660721134032),    // Central Indochina dry forests
        (142, 1.6920615104426),     // Central Korean deciduous forests
        (143, 2.16228419158302),    // Central Mexican matorral
        (144, 3.71613379296748),    // Central Pacific Northwest coastal forests
        (145, 0.989319566583888),   // Central Persian desert basins
        (147, 4.67883806623569),    // Central Range Papuan montane rain forests
        (148, 4.36706771936929),    // Papuan Central Range sub-alpine grasslands
        (149, 1.12328352330179),    // Central Ranges xeric scrub
        (151, 3.32096313919286),    // Central Tallgrass prairie
        (152, 0.826303795304645),   // Central Tibetan Plateau alpine steppe
        (153, 1.85971837254442),    // Central Zambezian wet miombo woodlands
        (154, 2.99044547132678),    // Cerrado
        (155, 2.02109875518812),    // Changbai Mountains mixed forests
        (156, 1.0856109607183),     // Changjiang Plain evergreen forests
        (157, 1.93668370741723),    // Chao Phraya freshwater swamp forests
        (158, 2.23920290887509),    // Chao Phraya lowland moist deciduous forests
        (160, -0.0210405194996651), // Cherskii-Kolyma mountain tundra
        (161, 0.679399978480928),   // Chhota-Nagpur dry deciduous forests
        (162, 2.60365692463316),    // Chiapas Depression dry forests
        (163, 2.45026587921876),    // Chiapas montane forests
        (164, 1.74110276918448),    // Chihuahuan desert
        (165, 1.44533707162547),    // Chilean Matorral
        (166, 2.86499943958747),    // Chimalapas montane forests
        (167, 1.30626124431965),    // Chin Hills-Arakan Yoma montane forests
        (168, 3.37921669850661),    // Chiquitano dry forests
        (169, 1.96157274748097),    // Chocó-Darién moist forests
        (171, -0.226029948351825),  // Chukchi Peninsula tundra
        (174, 1.68628890999041),    // Colorado Plateau shrublands
        (175, 2.38558024056),       // Colorado Rockies forests
        (176, 0.612106946218661),   // Comoros forests
        (177, 1.54534974024188),    // Cook Inlet taiga
        (178, 3.92816144889075),    // Cook Islands tropical moist forests
        (179, 1.28993835672988),    // Coolgardie woodlands
        (180, 0.799607956972179),   // Copper Plateau taiga
        (181, 2.92761152536065),    // Cordillera Central páramo
        (182, 4.44662637624143),    // Cordillera de Merida páramo
        (183, 3.44539446178338),    // Cordillera La Costa montane forests
        (184, 3.11460147696),       // Cordillera Oriental montane forests
        (185, 2.13384034844595),    // Corsican montane broadleaf and mixed forests
        (186, 2.38221717076262),    // Costa Rican seasonal moist forests
        (187, 1.65709588522591),    // Crete Mediterranean forests
        (188, 2.51853417753028),    // Crimean Submediterranean forest complex
        (189, -0.679891888578604),  // Cross-Niger transition forests
        (190, -0.684812199235456),  // Cross-Sanaga-Bioko coastal forests
        (191, 2.99006844776448),    // Cross-Timbers savanna-woodland
        (192, 2.6095354931862),     // Cuban cactus scrub
        (193, 2.59224388323921),    // Cuban dry forests
        (194, 2.63781391518964),    // Cuban moist forests
        (195, 2.62417320317419),    // Cuban pine forests
        (196, 2.51046372061792),    // Cuban wetlands
        (197, 1.92052163388829),    // Cyprus Mediterranean forests
        (198, 1.42695959857252),    // Da Hinggan-Dzhagdy Mountains conifer forests
        (199, 1.65596781802998),    // Daba Mountains evergreen forests
        (200, 1.49534243095636),    // Daurian forest steppe
        (201, -0.81310923982328),   // Davis Highlands tundra
        (202, 1.21385724353966),    // Deccan thorn scrub forests
        (203, 2.23979590776105),    // Dinaric Mountains mixed forests
        (204, 1.67200759885862),    // Djibouti xeric shrublands
        (205, 2.77220300724164),    // Drakensberg Escarpment savanna and thicket
        (206, 2.48360117490492),    // Drakensberg grasslands
        (207, -0.902197418475154),  // Rock and Ice
        (208, 2.45759638305879),    // Dry Chaco
        (209, 1.4925333169149),     // East Afghan montane conifer forests
        (210, 2.89134199476468),    // East African halophytics
        (211, 2.02036703422543),    // East African mangroves
        (212, 1.80859552499984),    // East African montane forests
        (214, 0.629875248094368),   // East Arabian fog shrublands and sand desert
        (215, 3.03364201593725),    // East Central Texas forests
        (216, 0.990575776775266),   // East Deccan dry-evergreen forests
        (217, 2.30727879315368),    // East European forest steppe
        (218, -0.128969809227929),  // East Sahara Desert
        (219, 0.877373864005263),   // East Saharan montane xeric woodlands
        (220, 0.666715629874607),   // East Siberian taiga
        (221, 1.0046661493695),     // East Sudanian savanna
        (222, 1.9330225366924),     // Eastern Anatolian deciduous forests
        (223, 2.01850136746443),    // Eastern Anatolian montane steppe
        (224, 2.1976955502959),     // Eastern Arc forests
        (225, 1.59056145548101),    // Eastern Australia mulga shrublands
        (226, 3.37718166383721),    // Eastern Australian temperate forests
        (227, 1.45406702972242),    // Eastern Canadian forests
        (228, 0.481407922553437),   // Eastern Canadian Shield taiga
        (229, 2.76336907871697),    // Eastern Cascades forests
        (230, 0.25757602010663),    // Eastern Congolian swamp forests
        (231, 2.33303573762273),    // Eastern Cordillera Real montane forests
        (232, 2.32091820009041),    // Eastern Canadian Forest-Boreal transition
        (233, 0.594664645324998),   // Eastern Gobi desert steppe
        (234, 3.13743115576298),    // Eastern Great Lakes lowland forests
        (235, 0.378574161626282),   // Eastern Guinean forests
        (236, 0.701791238332607),   // East Deccan moist deciduous forests
        (237, 2.05348134496644),    // Eastern Himalayan alpine shrub and meadows
        (238, 1.43432454034911),    // Eastern Himalayan broadleaf forests
        (239, 2.21959009724472),    // Eastern Himalayan subalpine conifer forests
        (240, 1.47547058724767),    // Eastern Java-Bali montane rain forests
        (241, 1.55963169704234),    // Eastern Java-Bali rain forests
        (242, 1.70330155604469),    // Eastern Mediterranean conifer-broadleaf forests
        (244, 2.36871800756633),    // Eastern Panamanian montane forests
        (245, 2.31591487830433),    // Ecuadorian dry forests
        (246, 2.29124520816482),    // Edwards Plateau savanna
        (247, 2.75415949437129),    // Einasleigh upland savanna
        (248, 1.99569474539137),    // Elburz Range forest steppe
        (251, 1.57247160538003),    // Emin Valley steppe
        (253, 2.26104841353324),    // English Lowlands beech forests
        (254, 2.55100634837897),    // Enriquillo wetlands
        (255, 1.44408973373155),    // Eritrean coastal desert
        (256, 1.75771804082622),    // Esperance mallee
        (257, 2.29638320851572),    // Espinal
        (258, 1.13408229841407),    // Ethiopian montane forests
        (259, 1.53138627234935),    // Ethiopian montane grasslands and woodlands
        (260, 1.49019870518834),    // Ethiopian montane moorlands
        (261, 2.20176921009748),    // Etosha Pan halophytics
        (262, 2.29890799752793),    // Euxine-Colchic broadleaf forests
        (263, 2.85302403565246),    // Everglades flooded grasslands
        (264, 1.71411435851463),    // Eyre and York mallee
        (265, 2.06799939768415),    // Faroe Islands boreal grasslands
        (267, 3.31598845365444),    // Fiji tropical dry forests
        (268, 3.44954705702171),    // Fiji tropical moist forests
        (269, 2.4059916823668),     // Fiordland temperate forests
        (270, 1.52367546810723),    // Flinders-Lofty montane woodlands
        (271, 3.33409014623309),    // Flint Hills tallgrass prairie
        (272, 2.11401295315825),    // Fraser Plateau and Basin conifer forests
        (273, 3.29461396066377),    // Fynbos shrubland
        (274, 1.63509389519312),    // Galápagos Islands xeric scrub
        (275, 1.38870035405295),    // Gariep Karoo
        (276, 1.11544449722041),    // Ghorat-Hazarajat alpine meadow
        (277, 0.728158449390612),   // Gibson desert
        (278, 2.07941351021723),    // Gissaro-Alai open woodlands
        (279, 0.834363688318575),   // Godavari-Krishna mangroves
        (280, 0.650377102794068),   // Gobi Lakes Valley desert steppe
        (282, 1.91633666285708),    // Great Basin montane forests
        (283, 1.44568265655722),    // Great Basin shrub steppe
        (284, 1.13433662077039),    // Great Lakes Basin desert steppe
        (285, 1.13014029481496),    // Great Sandy-Tanami desert
        (286, 0.710189022748993),   // Great Victoria desert
        (287, 2.939499600401),      // Greater Negros-Panay rain forests
        (288, 2.99095348726326),    // Guajira-Barranquilla xeric scrub
        (289, 2.29865627954794),    // Guayaquil flooded grasslands
        (290, 3.07581024233063),    // Guianan freshwater swamp forests
        (291, 3.11864921442801),    // Guianan Highlands moist forests
        (292, 2.99227586295268),    // Guianan lowland moist forests
        (293, 3.06062707419635),    // Guianan piedmont moist forests
        (294, 3.3754710265207),     // Guianan savanna
        (295, 0.123223148655336),   // Guinean forest-savanna
        (296, -0.705391707018958),  // Guinean mangroves
        (297, -0.481698430558482),  // Guinean montane forests
        (298, 1.50875720153397),    // Guizhou Plateau broadleaf and mixed forests
        (299, 0.87912736740615),    // Gulf of California xeric scrub
        (300, 2.74736164843066),    // Gulf of St. Lawrence lowland forests
        (301, 3.92215535076235),    // Gurupa várzea
        (302, 1.49291707583738),    // Hainan Island monsoon rain forests
        (303, 3.26314708845454),    // Halmahera rain forests
        (304, 1.09703405693258),    // Hampton mallee and woodlands
        (305, 3.31335758175202),    // Hawai'i tropical dry forests
        (306, 2.8792305442965),     // Hawai'i tropical high shrublands
        (307, 2.42973188791834),    // Hawai'i tropical low shrublands
        (308, 4.21950274892264),    // Hawai'i tropical moist forests
        (309, 1.13604671217165),    // Helanshan montane conifer forests
        (310, 2.50451452597255),    // Hengduan Mountains subalpine conifer forests
        (311, -1.00043675655423),   // Canadian High Arctic tundra
        (312, 1.23771181628564),    // High Monte
        (313, 2.30729406790799),    // Highveld grasslands
        (314, 0.872984167039178),   // Himalayan subtropical broadleaf forests
        (315, 1.72701572685116),    // Himalayan subtropical pine forests
        (316, 1.78540971142851),    // Hindu Kush alpine meadow
        (317, 2.543548651988),      // Hispaniolan dry forests
        (318, 2.59549860783145),    // Hispaniolan moist forests
        (319, 2.61080676748663),    // Hispaniolan pine forests
        (320, 1.22422059664122),    // Hobyo grasslands and shrublands
        (321, 2.62209194440329),    // Hokkaido deciduous forests
        (322, 2.70009492253786),    // Hokkaido montane conifer forests
        (323, 2.78593003046806),    // Honshu alpine conifer forests
        (324, 1.25558363317764),    // Horn of Africa xeric bushlands
        (325, 1.37079633280907),    // Huang He Plain mixed forests
        (326, 2.66991104548518),    // Humid Chaco
        (327, 2.62606551721736),    // Humid Pampas
        (328, 4.51387416513565),    // Huon Peninsula montane rain forests
        (329, 2.15686036138921),    // Iberian conifer forests
        (330, 2.1162063219789),     // Iberian sclerophyllous and semi-deciduous forests
        (331, 1.15529317775373),    // Iceland boreal birch forests and alpine tundra
        (332, 2.38650925748489),    // Illyrian deciduous forests
        (333, 1.90124585464254),    // Indochina mangroves
        (334, 1.16705718253232),    // Indus River Delta-Arabian Sea mangroves
        (335, 1.5056354672615),     // Indus Valley desert
        (336, 1.63509525399352),    // Inner Niger Delta flooded savanna
        (337, 0.542654059805047),   // Interior Alaska-Yukon lowland taiga
        (338, 3.71623359853051),    // Interior Plateau US Hardwood Forests
        (339, 0.348237039521708),   // Interior Yukon-Alaska alpine tundra
        (340, 2.97283857977522),    // Iquitos várzea
        (341, 2.19447420276354),    // Irrawaddy dry forests
        (342, 1.69825198211517),    // Irrawaddy freshwater swamp forests
        (343, 1.85473951726016),    // Irrawaddy moist deciduous forests
        (345, 2.34376534399956),    // Isthmian-Atlantic moist forests
        (346, 2.0268724379999),     // Isthmian-Pacific moist forests
        (347, 2.09051880454795),    // Italian sclerophyllous and semi-deciduous forests
        (348, 2.13235955066123),    // Itigi-Sumbu thicket
        (349, 3.00475707864291),    // Jalisco dry forests
        (350, 2.56354922171218),    // Jamaican dry forests
        (351, 2.6066308967211),     // Jamaican moist forests
        (352, 3.30536920599501),    // Japurá-Solimões-Negro moist forests
        (353, 2.64563996753698),    // Jarrah-Karri forest and shrublands
        (354, 1.16932664152046),    // Jian Nan subtropical evergreen forests
        (355, 0.7260445181391),     // Jos Plateau forest-grassland
        (357, 0.872136908693758),   // Junggar Basin semi-desert
        (358, 3.53198201749912),    // Juruá-Purus moist forests
        (359, -1.0537109845364),    // Kalaallit Nunaat High Arctic tundra
        (360, 2.08407506058852),    // Kalahari Acacia woodlands
        (361, 1.71857003766649),    // Kalahari xeric savanna
        (362, 1.13932698643581),    // Kamchatka-Kurile meadows and sparse forests
        (363, 1.21527956335479),    // Kamchatka taiga
        (364, 1.09611992758588),    // Kaokoveld desert
        (365, 1.28058463676695),    // Karakoram-West Tibetan Plateau alpine steppe
        (366, 2.49193221771175),    // Kayah-Karen montane rain forests
        (367, 1.94148156046138),    // Kazakh forest steppe
        (368, 0.950687256508894),   // Kazakh semi-desert
        (369, 1.56660422846225),    // Kazakh steppe
        (370, 1.42488803744517),    // Kazakh upland steppe
        (372, 1.14089858640893),    // Khangai Mountains alpine meadow
        (373, 1.27202359620705),    // Khangai Mountains conifer forests
        (374, 1.17761946727425),    // Khathiar-Gir dry deciduous forests
        (375, 2.17254877421709),    // Kimberly tropical savanna
        (376, 3.46061830918334),    // Kinabalu montane alpine meadows
        (377, 3.98735213178529),    // Klamath-Siskiyou forests
        (378, 3.61797839581279),    // Knysna-Amatole montane forests
        (379, 0.987960701755607),   // Kola Peninsula tundra
        (380, 1.27819539730574),    // Kopet Dag semi-desert
        (381, 1.40245704893428),    // Kopet Dag woodlands and forest steppe
        (382, 0.936566704780343),   // Kuh Rud and Eastern Iran montane woodlands
        (383, 2.83789353885811),    // Kwazulu Natal-Cape coastal forests
        (384, 2.86813725312923),    // La Costa xeric shrublands
        (385, 1.46136014810273),    // Lake Chad flooded savanna
        (386, 3.15488044610939),    // Lara-Falcón dry forests
        (387, 2.19483445154566),    // Leeward Islands moist forests
        (388, 1.84995900389168),    // Lesser Sundas deciduous forests
        (389, 2.71377707096684),    // Limpopo lowveld
        (390, 2.55377747449281),    // Llanos
        (391, 3.38318239872042),    // Lord Howe Island subtropical forests
        (392, 3.76567339878371),    // Louisiade Archipelago rain forests
        (393, 1.09568234879112),    // Low Monte
        (394, 0.410981666061132),   // Lower Gangetic Plains moist deciduous forests
        (395, 1.80749880337425),    // Luang Prabang montane rain forests
        (396, 2.26140966442979),    // Luzon montane rain forests
        (397, 2.39113431203096),    // Luzon rain forests
        (398, 2.28653128089192),    // Luzon tropical pine forests
        (399, 0.737763299506385),   // Madagascar dry deciduous forests
        (400, 0.807990787394825),   // Madagascar ericoid thickets
        (401, 0.574586924573433),   // Madagascar humid forests
        (402, 0.869316235697021),   // Madagascar mangroves
        (403, 1.11626906669234),    // Madagascar spiny thickets
        (404, 0.780380185333451),   // Madagascar subhumid forests
        (405, 0.912152900816372),   // Madagascar succulent woodlands
        (406, 3.81709435850547),    // Madeira-Tapajós moist forests
        (407, 0.848243433889278),   // Madeira evergreen forests
        (408, 2.49328197555643),    // Magdalena-Urabá moist forests
        (409, 2.72920116784395),    // Magdalena Valley dry forests
        (410, 2.94597801801554),    // Magdalena Valley montane forests
        (411, 2.09794223444819),    // Magellanic subpolar forests
        (412, 2.08030866232004),    // Makgadikgadi halophytics
        (413, 0.0345087559560982),  // Malabar Coast moist forests
        (414, 2.35908046281148), // Maldives-Lakshadweep-Chagos Archipelago tropical moist forests
        (416, 1.86087431583767), // Manchurian mixed forests
        (417, 1.12126507008311), // Mandara Plateau woodlands
        (418, 2.7865528793009),  // Maputaland coastal forests and woodlands
        (419, 3.11466814578366), // Maracaibo dry forests
        (420, 3.41851212938586), // Marajó várzea
        (421, 3.0755249305572),  // Maranhão Babaçu forests
        (422, 2.65412679291173), // Marañón dry forests
        (423, 3.33062842708322), // Marianas tropical dry forests
        (425, 2.18150601054913), // Masai xeric grasslands and shrublands
        (426, 0.842076652014313), // Mascarene forests
        (427, 3.75180865789535), // Mato Grosso tropical dry forests
        (428, 1.34626367129825), // Mediterranean Acacia-Argania dry woodlands and succulent thickets
        (429, 2.08138303785775), // Mediterranean conifer and mixed forests
        (430, 1.18462788328511), // Mediterranean dry woodlands and steppe
        (431, 1.74038286798916), // Mediterranean High Atlas juniper steppe
        (432, 1.7451607114598),  // Mediterranean woodlands and forests
        (433, 0.62022574392002), // Meghalaya subtropical forests
        (434, 2.23420703951312), // Mentawai Islands rain forests
        (435, 1.80259606933864), // Meseta Central matorral
        (436, 2.69522669730485), // Mesoamerican Gulf-Caribbean mangroves
        (437, 0.773517495760907), // Mesopotamian shrub desert
        (438, 1.75737488915731), // Mid-Canada Boreal Plains forests
        (439, -0.741850345105883), // Canadian Middle Arctic Tundra
        (440, 3.63857492504378), // Mid-Atlantic US coastal savannas
        (441, 1.27424183669018), // Midwest Canadian Shield forests
        (442, 3.32640425529744), // Mindanao-Eastern Visayas rain forests
        (443, 3.46054861504021), // Mindanao montane rain forests
        (444, 2.68259570371823), // Mindoro rain forests
        (445, 2.82316026106497), // Miskito pine forests
        (446, 3.6833998018299),  // Mississippi lowland forests
        (447, 1.6581094631319),  // Mitchell Grass Downs
        (448, 0.877132723438463), // Mizoram-Manipur-Kachin rain forests
        (449, 0.92357345921387), // Mojave desert
        (450, 1.19853780779889), // Mongolian-Manchurian grassland
        (451, 2.23495246121802), // Montana Valley and Foothill grasslands
        (452, 3.8572819394988),  // Monte Alegre várzea
        (453, 3.19827366666661), // Motagua Valley thornscrub
        (454, 0.236029062477386), // Mount Cameroon and Bioko montane forests
        (456, 1.76305666980133), // Murray-Darling woodlands and mallee
        (457, 1.08719441623081), // Muskwa-Slave Lake taiga
        (458, 1.83395050500403), // Myanmar Coast mangroves
        (459, 1.44297744289453), // Myanmar coastal rain forests
        (460, 1.86844956362743), // Nama Karoo shrublands
        (461, 1.29946249696077), // Namaqualand-Richtersveld steppe
        (462, 0.903733953654569), // Namib Desert
        (463, 1.5343111960327),  // Namibian savanna woodlands
        (464, 1.4404334544251),  // Nansei Islands subtropical evergreen forests
        (465, 1.83212166324668), // Napo moist forests
        (466, 2.52860886864489), // Naracoorte woodlands
        (467, 0.866706143877232), // Narmada Valley dry deciduous forests
        (468, 2.3865158640621),  // Nebraska Sand Hills mixed grasslands
        (469, 2.97396173761317), // Negro-Branco moist forests
        (470, 3.0511824447637),  // Nelson Coast temperate forests
        (471, 1.64501498153214), // Nenjiang River grassland
        (472, 3.93490632033621), // New Britain-New Ireland lowland rain forests
        (473, 3.8612438181251),  // New Britain-New Ireland montane rain forests
        (474, 2.31702610205349), // New Caledonia dry forests
        (475, 2.57461672453442), // New Caledonia rain forests
        (476, 2.66503326399862), // New England-Acadian forests
        (477, 4.44602431660313), // New Guinea mangroves
        (478, 3.20065974171009), // Nicobar Islands rain forests
        (479, -1.06047441968581), // Niger Delta swamp forests
        (480, 0.0140609761682934), // Nigerian lowland forests
        (481, 2.46673564118319), // Nihonkai evergreen forests
        (482, 2.78628644821277), // Nihonkai montane deciduous forests
        (483, 1.12482936678958), // Nile Delta flooded savanna
        (487, 0.518033679814302), // North Arabian desert
        (488, 0.444737563872754), // North Arabian highland shrublands
        (489, 3.51062829040445), // New Zealand North Island temperate forests
        (490, 0.259976423301978), // North Saharan Xeric Steppe and Woodland
        (491, 1.22442257184313), // Somali montane xeric woodlands
        (492, 0.60186092481888), // North Tibetan Plateau-Kunlun Mountains alpine desert
        (494, 0.397892631701345), // North Western Ghats moist deciduous forests
        (495, 0.0911511424767569), // North Western Ghats montane rain forests
        (496, 1.78408646979354), // Northeast China Plain deciduous forests
        (497, 1.27098961205068), // Northeast India-Myanmar pine forests
        (498, -0.352737753954282), // Northeast Siberian coastal tundra
        (499, 0.0946368918451431), // Northeast Siberian taiga
        (500, 3.03013285556584), // Northeast Brazil restingas
        (501, 0.30002470711584), // Northeast Congolian lowland forests
        (502, 2.40520394134289), // Northeast Himalayan subalpine conifer forests
        (503, 2.45479056228562), // Northeast Spain and Southern France Mediterranean forests
        (504, 3.60531502602429), // Northeast US Coastal forests
        (505, 2.28994823555976), // Northern Acacia-Commiphora bushlands and thickets
        (506, 2.06899635583688), // Northern Anatolian conifer and deciduous forests
        (507, 3.69732280860314), // Northern Andean páramo
        (508, 1.23040413088683), // Northern Annamites rain forests
        (509, 4.25009280634753), // Northern California coastal forests
        (510, 0.158199068148868), // Northern Canadian Shield taiga
        (511, 0.348709787481562), // Northern Congolian Forest-Savanna
        (512, 0.597788046400906), // North Deccan dry deciduous forests
        (513, 1.91424742580106), // Northern Indochina subtropical forests
        (514, 1.01441595368742), // Northern Khorat Plateau moist deciduous forests
        (515, 2.44420449602354), // Northern Mesoamerican Pacific mangroves
        (516, 4.60683808423761), // Northern New Guinea lowland rain and freshwater swamp forests
        (517, 4.59748302043145), // Northern New Guinea montane rain forests
        (518, 1.84097310363203), // Northern Shortgrass prairie
        (519, 1.94784248797897), // Northern Swahili coastal forests
        (520, 2.34078490339325), // Northern Tallgrass prairie
        (521, 2.20477355584549), // Northern Thailand-Laos moist deciduous forests
        (522, 1.31362865501416), // Northern Triangle subtropical forests
        (523, 1.82463799137378), // Northern Triangle temperate forests
        (524, 1.30824611626359), // Northern Vietnam lowland rain forests
        (525, 4.09859812476945), // Northland temperate kauri forests
        (526, 2.44518647139136), // Northwest Iberian montane forests
        (527, 0.32080278050485), // Northwest Russian-Novaya Zemlya tundra
        (528, 2.44677374377257), // Northwest Andean montane forests
        (529, 0.329290083840986), // Northwest Congolian lowland forests
        (531, 2.10892735029065), // Northwestern Himalayan alpine shrub and meadows
        (532, 1.4214907053853),  // Aravalli west thorn scrub forests
        (533, -0.698870787818633), // Novosibirsk Islands Arctic desert
        (534, 2.36249101708889), // Nujiang Langcang Gorge alpine conifer and mixed forests
        (535, 0.668550932609883), // Nullarbor Plains xeric shrublands
        (536, 2.61916785704806), // Nyanga-Chimanimani Montane forest-grassland
        (537, 2.84894237881618), // Oaxacan montane forests
        (539, 0.279322689433837), // Ogilvie-MacKenzie alpine tundra
        (540, 2.52358964429679), // Okanogan dry forests
        (541, 1.45626870600209), // Okhotsk-Manchurian taiga
        (542, 1.18029993246589), // Ordos Plateau steppe
        (543, 2.34423873808798), // Orinoco Delta swamp forests
        (544, 2.5902751172519),  // Orinoco wetlands
        (545, 0.80250123026353), // Orissa semi-evergreen forests
        (546, 3.7135203794168),  // Ozark Highlands mixed forests
        (547, 3.78088593761121), // Ozark Mountain forests
        (548, 1.38373487740093), // Pacific Coastal Mountain icefields and tundra
        (549, 4.06427710137302), // Palau tropical moist forests
        (550, 2.80254553378806), // Palawan rain forests
        (551, 1.95775652946116), // Palouse prairie
        (552, 1.38383961420843), // Pamir alpine desert and tundra
        (553, 2.5568139924928),  // Panamanian dry forests
        (554, 2.10987090413731), // Pannonian mixed forests
        (555, 3.21577715970136), // Pantanal
        (556, 2.50613937078649), // Pantanos de Centla
        (557, 3.20680933246307), // Pantepui forests & shrublands
        (558, 3.24240501815313), // Paraguaná xeric scrub
        (559, 2.78046268816694), // Paraná flooded savanna
        (560, 1.45657110791503), // Paropamisus xeric woodlands
        (561, 0.842408545711535), // Patagonian steppe
        (562, 2.07141483223078), // Patía valley dry forests
        (563, 2.90993925639631), // Peninsular Malaysian montane rain forests
        (564, 2.64833165883886), // Peninsular Malaysian peat swamp forests
        (565, 2.73238667359287), // Peninsular Malaysian rain forests
        (566, 2.40330043758462), // Pernambuco coastal forests
        (567, 2.55965396719744), // Pernambuco interior forests
        (568, 3.31283765403911), // Peruvian Yungas
        (569, 2.54967003853691), // Petén-Veracruz moist forests
        (570, 1.25066452358336), // Pilbara shrublands
        (571, 2.35275897149411), // Pindus Mountains mixed forests
        (572, 3.49349559413049), // Piney Woods
        (573, 2.45365351331289), // Po Basin mixed forests
        (574, 1.93473546439127), // Pontic steppe
        (576, 2.40982685455241), // Puerto Rican dry forests
        (577, 2.51631724609294), // Puerto Rican moist forests
        (578, 3.92550591680411), // Puget lowland forests
        (579, 3.81694584987305), // Purus-Madeira moist forests
        (580, 3.37570134513823), // Purus várzea
        (581, 2.7317854876217),  // Pyrenees conifer and mixed forests
        (582, 0.545156588159152), // Qaidam Basin semi-desert
        (583, 1.81830473353126), // Qilian Mountains conifer forests
        (584, 1.05080628504143), // Qilian Mountains subalpine meadows
        (585, 1.91444932855894), // Qin Ling Mountains deciduous forests
        (586, 2.62287594294425), // Qionglai-Minshan conifer forests
        (587, 3.78110526324193), // Queensland tropical rain forests
        (588, 2.68710428468354), // Rakiura Island temperate forests
        (589, 1.37124313168441), // Rann of Kutch seasonal salt marsh
        (591, 1.08635226191899), // Red River freshwater swamp forests
        (592, 0.71060782531543), // Red Sea-Arabian Desert shrublands
        (593, 0.768982976353252), // Red Sea coastal desert
        (594, 1.23505336102612), // Red Sea mangroves
        (595, 0.735321439887331), // Registan-North Pakistan sandy desert
        (596, 3.12220371734506), // Renosterveld shrubland
        (597, 3.01120321339419), // Richmond temperate forests
        (598, 3.25372857133801), // Rio Negro campinarana
        (599, 2.04182594978547), // Rodope montane mixed forests
        (600, 3.18633927977283), // Rwenzori-Virunga montane moorlands
        (601, 0.445515690257674), // Saharan halophytics
        (602, 1.12658518095111), // Sahelian Acacia savanna
        (603, 1.66786570316797), // Sakhalin Island taiga
        (604, 4.08694130273385), // Samoan tropical moist forests
        (606, 1.59223531954436), // San Lucan xeric scrub
        (607, 2.7962506046379),  // Santa Lucia Montane Chaparral & Woodlands
        (608, 3.36120995115351), // Santa Marta montane forests
        (609, 4.35871809460967), // Santa Marta páramo
        (610, -0.900933643617144), // São Tomé, Príncipe, and Annobón forests
        (611, 2.33474002654243), // Sarmatic mixed forests
        (612, 1.58834235855787), // Sayan alpine meadows and tundra
        (613, 1.73609733849323), // Sayan Intermontane steppe
        (614, 1.79105198124905), // Sayan montane conifer forests
        (615, 1.65535208371591), // Scandinavian and Russian taiga
        (616, 2.14827616402851), // Scandinavian coastal conifer forests
        (617, 1.26165253581043), // Scandinavian Montane Birch forest and grasslands
        (619, 1.49493194928192), // Sechura desert
        (620, 1.41467031890005), // Selenge-Orkhon forest steppe
        (621, 3.40258210633185), // Seram rain forests
        (623, 2.93362413125479), // Serra do Mar coastal forests
        (624, 1.59939520431131), // Sichuan Basin evergreen broadleaf forests
        (625, 1.76629113641035), // Sierra de la Laguna dry forests
        (626, 1.88547821048156), // Sierra de la Laguna pine-oak forests
        (627, 2.59010127129695), // Sierra de los Tuxtlas
        (628, 2.35006638573152), // Sierra Madre de Chiapas moist forests
        (629, 2.79610534569322), // Sierra Madre de Oaxaca pine-oak forests
        (630, 2.59370582612396), // Sierra Madre del Sur pine-oak forests
        (631, 2.77266492042951), // Sierra Madre Occidental pine-oak forests
        (632, 2.17935268948719), // Sierra Madre Oriental pine-oak forests
        (633, 3.06656816159283), // Sierra Nevada forests
        (634, 0.906382551578322), // Simpson desert
        (635, 3.05961046257355), // Sinaloan dry forests
        (636, 2.90693949869738), // Sinú Valley dry forests
        (637, 1.75527995367219), // Snake-Columbia shrub steppe
        (639, 0.760520541671438), // Socotra Island xeric shrublands
        (640, 2.70355546956888), // Solimões-Japurá moist forests
        (641, 3.64958596500681), // Solomon Islands rain forests
        (642, 1.56116729222719), // Somali Acacia-Commiphora bushlands and thickets
        (643, 2.3306126864249),  // Sonoran-Sinaloan subtropical dry forest
        (644, 1.23305230387712), // Sonoran desert
        (645, 2.05482642216341), // South American Pacific mangroves
        (647, 1.92756069939132), // South Apennine mixed montane forests
        (648, 2.3970275529548),  // South Central Rockies forests
        (649, 1.26994298882762), // South China-Vietnam subtropical evergreen forests
        (651, 1.47273607661752), // South Deccan Plateau dry deciduous forests
        (652, 1.00688668575983), // South Iran Nubo-Sindian desert and semi-desert
        (653, 2.30375089180584), // New Zealand South Island montane grasslands
        (654, 2.58717316399014), // New Zealand South Island temperate forests
        (656, 0.199491822838198), // South Sahara desert
        (657, 1.92651361457158), // South Siberian forest steppe
        (658, 1.66456754638787), // South Taiwan monsoon rain forests
        (660, 1.06352206886917), // South Western Ghats moist deciduous forests
        (661, 0.462486292363744), // South Western Ghats montane rain forests
        (662, 2.8157496551291),  // Southeast Australia temperate forests
        (663, 1.96333951292694), // Southeast Australia temperate savanna
        (664, 1.83856648625064), // Southeast Tibet shrublands and meadows
        (666, 3.31590316620417), // Southeast US conifer savannas
        (667, 2.00810411807053), // Southeast Iberian shrubs and woodlands
        (668, 1.58819890316697), // Southeast Indochina dry evergreen forests
        (669, 4.36195396360508), // Southeast Papuan rain forests
        (670, 2.24788569276994), // Southern Acacia-Commiphora bushlands and thickets
        (671, 2.80828401514624), // Southern Africa mangroves
        (672, 2.15400576531232), // Southern Anatolian montane conifer and deciduous forests
        (673, 1.15157198641438), // Southern Andean steppe
        (674, 2.94668331543804), // Southern Andean Yungas
        (675, 1.69483607764065), // Southern Annamites montane rain forests
        (676, 2.51591855849841), // Southern Atlantic Brazilian mangroves
        (677, 2.74046794211209), // Southern Cone Mesopotamian savanna
        (678, 0.749196394489163), // Southern Congolian forest-savanna
        (679, 3.49877001556552), // Southern Great Lakes forests
        (680, 0.925870085400845), // Southern Hudson Bay taiga
        (682, 1.46107725520158), // Southern Korea evergreen forests
        (683, 2.33928788885684), // Southern Mesoamerican Pacific mangroves
        (684, 4.41112398167224), // Southern New Guinea freshwater swamp forests
        (685, 4.60569148111836), // Southern New Guinea lowland rain forests
        (686, 2.71128631340105), // Southern Pacific dry forests
        (687, 1.84886191388646), // Southern Rift Montane forest-grassland
        (688, 2.03429850144532), // Southern Swahili coastal forests and woodlands
        (689, 2.01919507905306), // Southern Vietnam lowland dry forests
        (690, 3.53459788047351), // Southwest Amazon moist forests
        (691, 0.997559585089131), // Southwest Arabian Escarpment shrublands and woodlands
        (692, 0.986853518496514), // Southwest Arabian highland xeric scrub
        (693, 1.5601759705062),  // Southwest Australia savanna
        (694, 2.30023409446291), // Southwest Australia woodlands
        (695, 2.07036838242137), // Southwest Borneo freshwater swamp forests
        (696, 2.11611558132031), // Southwest Iberian Mediterranean sclerophyllous and mixed forests
        (697, 1.62097657342724), // Sri Lanka dry-zone dry evergreen forests
        (698, 1.32641738017634), // Sri Lanka lowland rain forests
        (699, 1.55303806474347), // Sri Lanka montane rain forests
        (701, 2.22126496849),    // Succulent Karoo xeric shrublands
        (702, 1.20917400556539), // Sudd flooded grasslands
        (703, 2.16687440330998), // Suiphun-Khanka meadows and forest meadows
        (704, 1.32125988902847), // Sulaiman Range alpine meadows
        (706, 2.42505855297097), // Sulawesi montane rain forests
        (707, 3.05068362613211), // Sulu Archipelago rain forests
        (708, 2.48486790000109), // Sumatran freshwater swamp forests
        (709, 2.54212659295363), // Sumatran lowland rain forests
        (710, 2.85709077365941), // Sumatran montane rain forests
        (711, 2.37630410067743), // Sumatran peat swamp forests
        (712, 3.1156420632812),  // Sumatran tropical pine forests
        (713, 1.83235702963103), // Sumba deciduous forests
        (714, 2.50452905616545), // Sunda Shelf mangroves
        (715, 2.21881186216311), // Sundaland heath forests
        (716, 0.383240366478809), // Sundarbans freshwater swamp forests
        (717, 0.217825493032596), // Sundarbans mangroves
        (718, 1.16242427979601), // Syrian xeric grasslands and shrublands
        (719, 2.17078387652739), // Taiheiyo evergreen forests
        (720, 2.68288438253134), // Taiheiyo montane deciduous forests
        (721, -0.392591293838957), // Taimyr-Central Siberian tundra
        (722, 1.90783821115595), // Taiwan subtropical evergreen forests
        (723, 0.303488831181414), // Taklimakan desert
        (724, 2.67755168709151), // Talamancan montane forests
        (725, 2.34081030263725), // Tamaulipan matorral
        (726, 1.98360425345205), // Tamaulipan mezquital
        (727, 3.87580727752368), // Tapajós-Xingu moist forests
        (728, 0.185783463410967), // Tarim Basin deciduous forests and steppe
        (729, 2.55776429011092), // Tasmanian Central Highland forests
        (730, 2.80764652822744), // Tasmanian temperate forests
        (731, 3.01442045227898), // Tasmanian temperate rain forests
        (732, 2.8720007413219),  // Tehuacán Valley matorral
        (733, 2.63065857355086), // Tenasserim-South Thailand semi-evergreen rain forests
        (734, 0.699005568683296), // Terai-Duar savanna and grasslands
        (735, 3.00966274518362), // Texas blackland prairies
        (736, 1.14360666932791), // Thar desert
        (737, 2.11067811058308), // Tian Shan foothill arid steppe
        (738, 1.85925904933655), // Tian Shan montane conifer forests
        (739, 1.60816371026826), // Tian Shan montane steppe and meadows
        (740, 0.0928248669157452), // Tibesti-Jebel Uweinat montane xeric woodlands
        (741, 1.08316648836198), // Tibetan Plateau alpine shrublands and meadows
        (742, 1.31488181927878), // Tigris-Euphrates alluvial salt marsh
        (743, 2.39831524133419), // Timor and Wetar deciduous forests
        (744, 0.904552744585701), // Tirari-Sturt stony desert
        (745, 3.4256046813779),  // Tongan tropical moist forests
        (746, 1.82725533273481), // Tonle Sap-Mekong peat swamp forests
        (747, 1.90028392190217), // Tonle Sap freshwater swamp forests
        (748, 0.0112991109222328), // Torngat Mountain tundra
        (749, 0.990045277658297), // Trans-Baikal Bald Mountain tundra
        (750, 1.60959357802196), // Trans-Baikal conifer forests
        (751, 2.8411273574324),  // Trans-Mexican Volcanic Belt pine-oak forests
        (752, 3.94005996330273), // Trans Fly savanna and grasslands
        (754, 2.14735985163171), // Trinidad and Tobago moist forest
        (756, 4.08217665340611), // Trobriand Islands rain forests
        (759, 2.06580879234965), // Tumbes-Piura dry forests
        (760, 1.93220302636292), // Tyrrhenian-Adriatic sclerophyllous and mixed forests
        (761, 3.46225854923671), // Uatumã-Trombetas moist forests
        (762, 3.00218272101145), // Ucayali moist forests
        (763, 1.05894205821123), // Upper Gangetic Plains moist deciduous forests
        (764, 3.06789680518265), // Upper Midwest US forest-savanna transition
        (765, 1.64620205337754), // Urals montane forest and taiga
        (766, 2.71382368973967), // Uruguayan savanna
        (767, 2.08552642094198), // Ussuri broadleaf and mixed forests
        (768, 2.81909946711931), // Valdivian temperate forests
        (769, 3.03412273704107), // Vanuatu rain forests
        (770, 3.76829368618187), // Venezuelan Andes montane forests
        (771, 2.98965079137992), // Veracruz dry forests
        (772, 2.87064883779158), // Veracruz moist forests
        (773, 2.91201180430084), // Veracruz montane forests
        (774, 1.91693078339478), // Victoria Plains tropical savanna
        (775, 4.13982969377666), // Vogelkop-Aru lowland rain forests
        (776, 4.2060228798854),  // Vogelkop montane rain forests
        (777, 2.41298985647148), // Wasatch and Uinta montane forests
        (778, 0.784741307331435), // Watson Highlands taiga
        (779, -0.035543277888458), // West Sahara desert
        (780, 0.0823301514042148), // West Saharan montane xeric woodlands
        (781, 1.24175460987456), // West Siberian taiga
        (782, 0.885068187281528), // West Sudanian savanna
        (783, 1.02065387271679), // Western Australian Mulga shrublands
        (784, 0.353993503352529), // Western Congolian swamp forests
        (785, 1.61905710945639), // Western Ecuador moist forests
        (786, 2.33691664250799), // Western European broadleaf forests
        (787, 2.65863589760327), // Western Great Lakes forests
        (788, -0.602810976454274), // Western Guinean lowland forests
        (789, 2.81953543855901), // Western Gulf coastal grasslands
        (790, 1.68204500589665), // Western Himalayan alpine shrub and meadows
        (791, 2.22295063382383), // Western Himalayan broadleaf forests
        (792, 2.0721609037739),  // Western Himalayan subalpine conifer forests
        (793, 1.84231150337042), // Western Java montane rain forests
        (794, 1.8297589849621),  // Western Java rain forests
        (796, 2.03811211895464), // Western shortgrass prairie
        (797, 1.97039257241087), // Western Siberian hemiboreal forests
        (798, 2.92610534359058), // Westland temperate forests
        (799, 4.36846550403486), // Willamette Valley oak savanna
        (800, 1.99240539198421), // Windward Islands moist forests
        (801, -0.485377927319135), // Wrangel Island Arctic desert
        (802, 1.61631446211846), // Wyoming Basin shrub steppe
        (803, 3.66730081576434), // Xingu-Tocantins-Araguaia moist forests
        (804, -0.0904575520791858), // Yamal-Gydan tundra
        (806, 4.42438550989839), // Yapen rain forests
        (807, 1.34259170030904), // Yarlung Zanbo arid steppe
        (808, 1.20877279680191), // Yellow Sea saline meadow
        (809, 2.47786084669482), // Yucatán dry forests
        (810, 2.63570573799733), // Yucatán moist forests
        (811, 2.11011336457414), // Yunnan Plateau subtropical evergreen forests
        (812, 1.60488881354278), // Zagros Mountains forest steppe
        (813, 2.4768724209414),  // Zambezian-Limpopo mixed woodlands
        (814, 2.32064208537758), // Zambezian Baikiaea woodlands
        (815, 2.12900002259769), // Zambezian coastal flooded savanna
        (816, 2.02709419159436), // Dry miombo woodlands
        (817, 2.17700353126079), // Zambezian evergreen dry forests
        (818, 2.07555149362847), // Zambezian flooded grasslands
        (819, 2.47794413711992), // Zambezian mopane woodlands
        (820, 2.7579139818189),  // Brazilian Atlantic dry forests
        (821, 1.18621978928392), // Northern Cordillera forests
        (823, 1.3989810819149),  // Victoria Basin forest-savanna
        (824, 0.72331629018183), // Western Congolian forest-savanna
        (825, 2.27121527834918), // East African montane moorlands
        (826, 3.03893208612349), // Queen Charlotte Islands conifer forests
        (827, 2.35387728613581), // Northern Pacific Alaskan coastal forests
        (828, 0.226082963697568), // Northwest Territories taiga
        (829, 0.981120166089835), // Southwest Arabian coastal xeric shrublands
        (830, 0.695035057538745), // Arabian-Persian Gulf coastal plain desert
        (831, 0.595007357444992), // South Arabian plains and plateau desert
        (832, 0.54491281744585), // Arabian desert
        (833, 2.51078284521512), // Bahamian pineyards
        (835, 0.862922606189755), // Kamchatka tundra
        (836, 3.49812738805455), // Tocantins/Pindare moist forests
        (837, -0.410918609811858), // Kalaallit Nunaat Arctic steppe
        (838, 3.05259615282773), // North Cascades conifer forests
        (839, 2.53959368368306), // Northern Rockies conifer forests
        (841, 1.09173690242668), // Southwest Arabian montane woodlands and grasslands
        (842, 0.852465759274842), // South Arabian fog woodlands, shrublands, and dune
        (843, 2.4379215346006),  // Trinidad and Tobago dry forest
        (844, 2.03442086986203), // Lesser Antillean dry forests
        (845, 2.59516415454447), // North Atlantic moist mixed forests
        (848, 2.50804814375686), // Sulawesi lowland rain forests
    ]
    .drain(..)
    .collect();

    match logpopdensity.get(&ecoregion) {
        None => return 0.,
        Some(ldensity) => return area_in_100_km2 * ldensity.exp(),
    }
}

#[derive(Clone, Copy)]
pub struct Ecovector {
    pub entries: [f32; ATTESTED_ECOREGIONS],
}

use std::ops::{Index, Mul};

impl Mul<f32> for Ecovector {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        let mut result = [0.; ATTESTED_ECOREGIONS];
        for (i, item) in self.entries.iter().enumerate() {
            result[i] = f32::max(0.5, item * rhs); //TODO: Move this part of the logic elsewhere
        }
        Ecovector { entries: result }
    }
}

impl Index<usize> for Ecovector {
    type Output = f32;

    fn index(&self, index: usize) -> &f32 {
        &self.entries[index]
    }
}

impl Default for Ecovector {
    fn default() -> Self {
        Ecovector {
            entries: [0.5; ATTESTED_ECOREGIONS],
        }
    }
}
