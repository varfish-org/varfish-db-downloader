# Changelog

## [1.0.0](https://www.github.com/bihealth/varfish-db-downloader/compare/v0.3.0...v1.0.0) (2024-01-19)


### ⚠ BREAKING CHANGES

* refactoring for bollonaster (#22)

### Features

* add DECIPHER HI predictions v3 ([#63](https://www.github.com/bihealth/varfish-db-downloader/issues/63)) ([375dd36](https://www.github.com/bihealth/varfish-db-downloader/commit/375dd362d53a61a8b0db13ef492f71766ab06638))
* add missing background data for GRCh38 ([#5](https://www.github.com/bihealth/varfish-db-downloader/issues/5)) ([#37](https://www.github.com/bihealth/varfish-db-downloader/issues/37)) ([9ce1cab](https://www.github.com/bihealth/varfish-db-downloader/commit/9ce1cab9e94b34f349bc1752c4cd7b3d1d712cfb))
* add nightly test for upstream URL availability ([#6](https://www.github.com/bihealth/varfish-db-downloader/issues/6)) ([#32](https://www.github.com/bihealth/varfish-db-downloader/issues/32)) ([e30b093](https://www.github.com/bihealth/varfish-db-downloader/commit/e30b093e9b2ac87960269028e1e57009247b2dab))
* add reduced development dataset ([#44](https://www.github.com/bihealth/varfish-db-downloader/issues/44)) ([#47](https://www.github.com/bihealth/varfish-db-downloader/issues/47)) ([dd8170e](https://www.github.com/bihealth/varfish-db-downloader/commit/dd8170e55df207e2258aef3e4796782d3c868020))
* adding AlphaMissense scores ([#60](https://www.github.com/bihealth/varfish-db-downloader/issues/60)) ([53e4209](https://www.github.com/bihealth/varfish-db-downloader/commit/53e4209e79f967c5131f05b3914f37c30f086819))
* adding annonars functional ([#68](https://www.github.com/bihealth/varfish-db-downloader/issues/68)) ([f00db04](https://www.github.com/bihealth/varfish-db-downloader/commit/f00db04286b6406346434be2fee4485740f84a22))
* adding clingen curation data ([#55](https://www.github.com/bihealth/varfish-db-downloader/issues/55)) ([51435df](https://www.github.com/bihealth/varfish-db-downloader/commit/51435dfd86bc4f02ae4ae71cb97ae3a29d54e695))
* adding DECIPHER HI to annonars genes ([#69](https://www.github.com/bihealth/varfish-db-downloader/issues/69)) ([307cf58](https://www.github.com/bihealth/varfish-db-downloader/commit/307cf5838d56d864f733cce27400269f4866848b))
* adding files for varfish-server-worker ([#50](https://www.github.com/bihealth/varfish-db-downloader/issues/50)) ([#51](https://www.github.com/bihealth/varfish-db-downloader/issues/51)) ([ba84443](https://www.github.com/bihealth/varfish-db-downloader/commit/ba84443ff611f21341414b4b87a55d55f7f95a45))
* adding pHaplo, pTriplo, sHet as seen in DECIPHER ([#56](https://www.github.com/bihealth/varfish-db-downloader/issues/56)) ([b0d6884](https://www.github.com/bihealth/varfish-db-downloader/commit/b0d6884c949899f154dd7029a9e59cabcaed68c4))
* adding test mode ([#7](https://www.github.com/bihealth/varfish-db-downloader/issues/7)) ([#30](https://www.github.com/bihealth/varfish-db-downloader/issues/30)) ([a415bfe](https://www.github.com/bihealth/varfish-db-downloader/commit/a415bfe64f8e8a21f76e9a2dc34180645244dbcf))
* adjusting output paths ([#42](https://www.github.com/bihealth/varfish-db-downloader/issues/42)) ([4b01fe9](https://www.github.com/bihealth/varfish-db-downloader/commit/4b01fe9340bfac15883fa217702997a4e10499e8))
* annotating with Orphanet diseases ([#58](https://www.github.com/bihealth/varfish-db-downloader/issues/58)) ([7d43865](https://www.github.com/bihealth/varfish-db-downloader/commit/7d43865a48cf68342a1f1f23ab33ffeee5a22a78))
* annotation of genes with OMIM diseases ([#57](https://www.github.com/bihealth/varfish-db-downloader/issues/57)) ([6c19f59](https://www.github.com/bihealth/varfish-db-downloader/commit/6c19f5927870868b43d651b8eb6356e527a494f5))
* binary conversion of sequence dbs with annonars and worker ([#35](https://www.github.com/bihealth/varfish-db-downloader/issues/35)) ([#40](https://www.github.com/bihealth/varfish-db-downloader/issues/40)) ([8d1542f](https://www.github.com/bihealth/varfish-db-downloader/commit/8d1542f0e4ff312f4e03828326ac7237f5c96f50))
* build annonars regions (clingen dosage) ([#67](https://www.github.com/bihealth/varfish-db-downloader/issues/67)) ([664295d](https://www.github.com/bihealth/varfish-db-downloader/commit/664295d8db03130bb0a3b52520f17b166874e40b))
* build output files for background data containers ([#35](https://www.github.com/bihealth/varfish-db-downloader/issues/35)) ([#41](https://www.github.com/bihealth/varfish-db-downloader/issues/41)) ([f3d936c](https://www.github.com/bihealth/varfish-db-downloader/commit/f3d936c2e970fed49b57b34d161174b2fc4b1684))
* bump dbvar to 20231030 ([#81](https://www.github.com/bihealth/varfish-db-downloader/issues/81)) ([4fe6866](https://www.github.com/bihealth/varfish-db-downloader/commit/4fe68668c27ef8b0796a5ddf18f0aa6b368b4ae8))
* cleanup and restructuring of output and code ([#34](https://www.github.com/bihealth/varfish-db-downloader/issues/34)) ([#36](https://www.github.com/bihealth/varfish-db-downloader/issues/36)) ([400119f](https://www.github.com/bihealth/varfish-db-downloader/commit/400119fc90dc48662a50ccf427f23e1e2e9f8d54))
* download of ClinGen dosage sensitivity ([#61](https://www.github.com/bihealth/varfish-db-downloader/issues/61)) ([17fa976](https://www.github.com/bihealth/varfish-db-downloader/commit/17fa976c4075b412cea40e1bb537fba0e56c8ce6))
* download of ClinGen dosage sensitivity ([#62](https://www.github.com/bihealth/varfish-db-downloader/issues/62)) ([b6bc6c5](https://www.github.com/bihealth/varfish-db-downloader/commit/b6bc6c56824dcbdda306f3cebeca39b4a74c263d))
* import of gnomAD SV data into RocksDB ([#66](https://www.github.com/bihealth/varfish-db-downloader/issues/66)) ([3fd72dd](https://www.github.com/bihealth/varfish-db-downloader/commit/3fd72ddd835d38159af18df36c4119105b061587))
* importing GTex data into annonars genes database ([#59](https://www.github.com/bihealth/varfish-db-downloader/issues/59)) ([a2e63f0](https://www.github.com/bihealth/varfish-db-downloader/commit/a2e63f0981268043c452abc0da1928be0e4f2128))
* integrating gnomAD v4 ([#75](https://www.github.com/bihealth/varfish-db-downloader/issues/75)) ([d529d7e](https://www.github.com/bihealth/varfish-db-downloader/commit/d529d7e33fc85c352924a06eb8e5d075099ee339))
* integrating PanelApp download for annonars ([#79](https://www.github.com/bihealth/varfish-db-downloader/issues/79)) ([#80](https://www.github.com/bihealth/varfish-db-downloader/issues/80)) ([3a332a0](https://www.github.com/bihealth/varfish-db-downloader/commit/3a332a078ffa9f3e9cb4de927e98ffa12496a5a4))
* introduce comprehensive spec YAML files ([#33](https://www.github.com/bihealth/varfish-db-downloader/issues/33)) ([#43](https://www.github.com/bihealth/varfish-db-downloader/issues/43)) ([3f70b57](https://www.github.com/bihealth/varfish-db-downloader/commit/3f70b5775cdbce77a28e187ead6ccec3dfadafa6))
* make runnable with Slurm ([#38](https://www.github.com/bihealth/varfish-db-downloader/issues/38)) ([#39](https://www.github.com/bihealth/varfish-db-downloader/issues/39)) ([2f6deac](https://www.github.com/bihealth/varfish-db-downloader/commit/2f6deacbe222873eeb7758ca08b662671552a948))
* provide enrichment probesets for Agilent/IDT/Twist ([#52](https://www.github.com/bihealth/varfish-db-downloader/issues/52)) ([#53](https://www.github.com/bihealth/varfish-db-downloader/issues/53)) ([c156f7f](https://www.github.com/bihealth/varfish-db-downloader/commit/c156f7f1d567a37e9eaf691e2b0a511e1f3b4d2b))
* refactoring for bollonaster ([#22](https://www.github.com/bihealth/varfish-db-downloader/issues/22)) ([1b35484](https://www.github.com/bihealth/varfish-db-downloader/commit/1b354843482bd4db72c45bb87d5d9143d25b3b87))
* replace orphapacket by orphadata API access ([#84](https://www.github.com/bihealth/varfish-db-downloader/issues/84)) ([#85](https://www.github.com/bihealth/varfish-db-downloader/issues/85)) ([482bdd0](https://www.github.com/bihealth/varfish-db-downloader/commit/482bdd0f4982172d1548d571a6b9f60e4da0ef29))
* retrieve downloader version via cli ([#83](https://www.github.com/bihealth/varfish-db-downloader/issues/83)) ([06c800c](https://www.github.com/bihealth/varfish-db-downloader/commit/06c800c8336a60723e2aad9f8b4e93fdecb38d7e))
* updating annonars genes with ClinGen dosage & DOMINO ([#65](https://www.github.com/bihealth/varfish-db-downloader/issues/65)) ([0797f2d](https://www.github.com/bihealth/varfish-db-downloader/commit/0797f2dcfd3e0d639c6c0648b65ac326df8c81d7))
* upgrade dbNSFP to v4.5 ([#77](https://www.github.com/bihealth/varfish-db-downloader/issues/77)) ([#78](https://www.github.com/bihealth/varfish-db-downloader/issues/78)) ([62d2f9a](https://www.github.com/bihealth/varfish-db-downloader/commit/62d2f9a2e2fe7a221c7cf1648bd88c01cbe656bc))
* version bump of mehari and transcripts ([#64](https://www.github.com/bihealth/varfish-db-downloader/issues/64)) ([bbcf2f2](https://www.github.com/bihealth/varfish-db-downloader/commit/bbcf2f2f685fa81e86091282185fac29aa1cd3f4))


### Bug Fixes

* adjusting to GTEx upstream URL change ([#71](https://www.github.com/bihealth/varfish-db-downloader/issues/71)) ([b3ec6d8](https://www.github.com/bihealth/varfish-db-downloader/commit/b3ec6d8da5a5eda339894644ae7efee6b95daede))
* removing duplicates in conditions JSONL ([#86](https://www.github.com/bihealth/varfish-db-downloader/issues/86)) ([16edd9b](https://www.github.com/bihealth/varfish-db-downloader/commit/16edd9b305bc2b16d9a79aa4686ff0f8cb3fb4a0))
* small adjustments to reduce data list ([#48](https://www.github.com/bihealth/varfish-db-downloader/issues/48)) ([e899051](https://www.github.com/bihealth/varfish-db-downloader/commit/e899051d01685b8b0c5720fb1973582b8f675742))
