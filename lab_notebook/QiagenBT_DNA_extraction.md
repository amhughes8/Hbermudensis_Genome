# Qiagen Blood and Tissue DNA extractions for Illumina + ONT
## Date: November 3, 2025
## Who: Annabel Hughes

1. Using gill and fin tissue from 5 blue angelfish collected in September 2025 from Texas. All tissues were preserved in either DNA/RNA Shield or RNAlater to see which gives us the best quality extraction:

| Tissue Type | Preservation Method | Freezer Storage | Sample ID |
|-----------|-------------------|----------------|----------|
| fin | RNAlater | -80C | HBE_TEX_092501_F |
| gill | RNAlater | -80C | HBE_TEX_092501_G |
| fin | DNA/RNA Shield | -80C | HBE_TEX_092502_F |
| gill | DNA/RNA Shield | -80C | HBE_TEX_092502_G |
| fin | DNA/RNA Shield | -80C | HBE_TEX_092503_F |
| gill | DNA/RNA Shield | -80C | HBE_TEX_092503_G |
| fin | RNAlater | -80C | HBE_TEX_092504_F |
| gill | RNAlater | -80C | HBE_TEX_092504_G |
| fin | DNA/RNA Shield | -80C | HBE_TEX_092505_F |
| gill | DNA/RNA Shield | -80C | HBE_TEX_092505_G |

2. Incubated for 2 hr 10 min at 56C (flicked and spun every 20-30 mins).
   - tissues stored in DNA/RNA shield were very slimy and took way longer to lyse. I probably could have let them go longer, because they ended up clogging the filters from the kit, but the tissues stored in RNAlater were lysed within 1 hr so I didn't want to leave them on the heat for too long.
3. Eluted in 200 uL EB (Tris-HCL pH=8.2). The EB was not pre-heated.
4. Repeated elution with original 200 uL

## Qubit Results:
| Sample ID | Concentration [ng/uL] |
|---------|---------------------|
| HBE_TEX_092501_F | 125 |
| HBE_TEX_092501_G | 38.8 |
| HBE_TEX_092502_F | 2.26 |
| HBE_TEX_092502_G | 19.4 |
| HBE_TEX_092503_F | 1.74 |
| HBE_TEX_092503_G | 2.54 |
| HBE_TEX_092504_F | 8.88 |
| HBE_TEX_092504_G | 41.2 |
| HBE_TEX_092505_F | 3.9 |
| HBE_TEX_092505_G | 1.31 |

## Nanodrop Results:
| Sample ID | Concentration [ng/uL] | 260/280 | 260/230 |
|---------|--------------------|---------|---------|
| HBE_TEX_092501_F | 183.7 | 1.98 | 2.25 |
| HBE_TEX_092501_G | 160 | 1.97	| 2.18 |
| HBE_TEX_092502_F | 27 | 2.24 | 1.41 |
| HBE_TEX_092502_G | 27| 1.81 |	1.11 |
| HBE_TEX_092503_F | 44 | 2.11 | 1.1 |
| HBE_TEX_092503_G | 15 | 1.96 | 0.19 |
| HBE_TEX_092504_F | 59.6 | 1.98 | 2.02 |
| HBE_TEX_092504_G | 238.5 | 1.91 |	2.26 |
| HBE_TEX_092505_F | 43.6 | 2.04 |	1.23 |
| HBE_TEX_092505_G | 8.8 | 1.76 |	0.36 |

## Tapestation Results:
Ran a gDNA tape with HBE_TEX_092501_F, HBE_TEX_092501_G, and HBE_TEX_092504_G. HBE_TEX_092501_F looked the best through all quality checks:

![plot](https://github.com/amhughes8/Hbermudensis_Genome/blob/main/photos/tapestation_results_HBE_gDNA.png)
