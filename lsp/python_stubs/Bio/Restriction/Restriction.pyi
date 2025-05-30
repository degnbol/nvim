from Bio.Restriction.PrintFormat import PrintFormat
from Bio.Seq import Seq
from _typeshed import Incomplete

__all__ = ['FormattedSeq', 'Analysis', 'RestrictionBatch', 'AllEnzymes', 'CommOnly', 'NonComm', 'PvuII', 'CcrNAIII', 'GauT27I', 'MspI7II', 'Sag901I', 'UbaF13I', 'SimI', 'BspD6I', 'WviI', 'BtuMI', 'BcnI', 'AgeI', 'KroI', 'BstSFI', 'HhaI', 'Hpy99I', 'PssI', 'AclWI', 'BsnI', 'AspJHL3II', 'Bbr52II', 'Dpi3090II', 'PgaP73III', 'Kor51II', 'Ssp714II', 'ApyPI', 'YkrI', 'PspPPI', 'Ple19I', 'OliI', 'HspAI', 'Bsp120I', 'SfuI', 'Tsp45I', 'DseDI', 'PshAI', 'Cco14983V', 'FspPK15I', 'Rtr1953I', 'Mlu211III', 'UbaF11I', 'BsiI', 'BcefI', 'TsoI', 'BstSNI', 'AclI', 'HpySE526I', 'AvaII', 'BstPI', 'FseI', 'DrdI', 'XhoII', 'Nli3877I', 'AarI', 'Bsh1236I', 'Jsp2502II', 'Acc65V', 'BanLI', 'Dpi3069I', 'SpnRII', 'PfrJS12V', 'Hin4II', 'MvnI', 'HpaII', 'Bsp19I', 'SaqAI', 'TfiI', 'Psp6I', 'PacI', 'CaiI', 'SsiI', 'NruI', 'Rsp531II', 'Cco11366VI', 'Fnu11326II', 'MkaDII', 'TsuI', 'AspBHI', 'TkoI', 'BstFNI', 'AccIII', 'Hin1I', 'BstDSI', 'BstXI', 'AspS9I', 'BmtI', 'BspMII', 'HgiJII', 'TaqII', 'Bse8I', 'AccIX', 'Awo1030IV', 'CspL61I', 'Sno506I', 'PflPt14I', 'Jma19592I', 'StsI', 'Eco57MI', 'Msp20I', 'HinP1I', 'BsiSI', 'StyI', 'KspI', 'PteI', 'PfoI', 'BstMWI', 'Bst2BI', 'NaeI', 'CcaP7V', 'Fco1691IV', 'Rsp008IV', 'Mch946II', 'PsuGI', 'SstE37I', 'FbaI', 'BssT1I', 'Bsp1286I', 'XmaI', 'AjnI', 'ApaI', 'TagI', 'FmuI', 'MnlI', 'BoxI', 'AteTI', 'Csp2014I', 'Hso63373III', 'SmaUMH8I', 'Pfl3756II', 'Bce83I', 'MroXI', 'BshTI', 'PspXI', 'FatI', 'SinI', 'Bst4CI', 'Hin1II', 'PfeI', 'BauI', 'PpsI', 'MscI', 'Rmu369III', 'Lsp6406VI', 'Cba16038I', 'Esp3007I', 'PenI', 'RpaB5I', 'XbaI', 'Eco52I', 'AccI', 'Bsp1720I', 'Bsc4I', 'SciI', 'ApaBI', 'PleI', 'MboII', 'TspDTI', 'BmcAI', 'Sma10259II', 'Pdi8503III', 'Asp114pII', 'Cre7908I', 'HpyUM032XIV', 'MluNI', 'EagI', 'BsePI', 'Psp1406I', 'ScrFI', 'NmuCI', 'EcoT22I', 'BsiEI', 'HpaI', 'MspJI', 'BssSI', 'Cau10061II', 'ErhG4T10I', 'Rho5650I', 'LsaDS4I', 'CjuII', 'RpaI', 'Hin4I', 'TspMI', 'CspAI', 'BshNI', 'BanII', 'NspBII', 'TspEI', 'UnbI', 'Asp700I', 'Lsp1109I', 'HphI', 'Mva1269I', 'AspNIH4III', 'Cpe10578V', 'HpyPU007XIX', 'SgrAII', 'PcaII', 'MalI', 'DpnII', 'Bse118I', 'PscI', 'RsrII', 'MspR9I', 'BstNSI', 'BseSI', 'BbvCI', 'LpnPI', 'HincII', 'RflFIII', 'Lpn12272I', 'Cal14237I', 'Ehi46392I', 'BspGI', 'RdeGBII', 'Bsp24I', 'Sse9I', 'CciNI', 'BseBI', 'AhdI', 'LpnI', 'SplI', 'SfeI', 'MbiI', 'FokI', 'BstF5I', 'GsuI', 'AfeI', 'Sfl13829III', 'Pba2294I', 'AspAMDIV', 'Cly7489II', 'HpyLIM6XII', 'KroNI', 'CviAII', 'Bpu14I', 'Pfl23II', 'PspGI', 'MabI', 'BstHHI', 'BmeRI', 'FspI', 'FspEI', 'PsrI', 'RdeGBI', 'Lpl1004II', 'BspNCI', 'Abr4036II', 'EcoMVII', 'AlwFI', 'PspPRI', 'SgrAI', 'BsuTUI', 'BpuMI', 'AhaIII', 'Ppu10I', 'SauI', 'CchII', 'BmgBI', 'Eco31I', 'BspCNI', 'BtsCI', 'AcvI', 'HpyAS001VI', 'Sep11964I', 'Asp103I', 'Cko11077IV', 'PaePA99III', 'Hpy166II', 'Cfr10I', 'BfrI', 'PauI', 'NciI', 'HpyF3I', 'BspOI', 'Bbv12I', 'FaqI', 'BarI', 'DraI', 'AbaCIII', 'Bsp460III', 'EcoE1140I', 'Ran11014IV', 'Lmo370I', 'CdiI', 'PlaDI', 'ZrmI', 'SalI', 'BstX2I', 'BmrFI', 'BetI', 'Hpy178III', 'ChaI', 'AccBSI', 'BtgZI', 'BseRI', 'BtsI', 'AanI', 'Aod1I', 'CjeNII', 'SenSARA26III', 'Pae8506I', 'HpyAXVI_mut2', 'GluI', 'FspAI', 'Bsu15I', 'HinfI', 'PagI', 'BcuI', 'Alw21I', 'BfoI', 'CseI', 'AjuI', 'Cac8I', 'Bps6700III', 'Aba13301I', 'Eco9699II', 'Pst273I', 'LlaG50I', 'Yps3606I', 'NlaCI', 'Zsp2I', 'StuI', 'PaeR7I', 'BstBI', 'BmeT110I', 'DsaI', 'BsmFI', 'BpuEI', 'BsrDI', 'HpyAXIV', 'AhyRBAHI', 'CjeFIII', 'Sen17963III', 'PacIII', 'NgoAVII', 'UpaP162I', 'EheI', 'BsrFI', 'AdeI', 'ErhI', 'Fnu4HI', 'NdeII', 'AsuNHI', 'BstV2I', 'BstUI', 'Vtu19109I', 'PspR84I', 'Ble402II', 'Eco9020I', 'LcrJM4II', 'DrdV', 'TstI', 'SrfI', 'NheI', 'BstAFI', 'TaiI', 'Bme18I', 'XmaJI', 'CauII', 'MlyI', 'BsmAI', 'BmrI', 'BseNI', 'AcoY31II', 'Cje265V', 'Sen5794III', 'OgrI', 'Hpy99XXII', 'EsaBC3I', 'EcoO65I', 'EcoO109I', 'MspCI', 'Eco53kI', 'BsiWI', 'AspA2I', 'AasI', 'McaTI', 'BstMAI', 'BsaAI', 'Bhe175II', 'AlfI', 'Eco8164I', 'Vdi96II', 'PspD7DII', 'Lcr047I', 'DraRI', 'PpiI', 'Sse8387I', 'BfmI', 'SmiI', 'NdeI', 'BssNI', 'VneI', 'SetI', 'Bpu10I', 'BseGI', 'Cgl13032I', 'Hpy99XIV', 'Sdy7136I', 'NspES21II', 'TscAI', 'SbfI', 'Eco147I', 'BsaWI', 'Eco91I', 'DdeI', 'MroI', 'AsiGI', 'HauII', 'FalI', 'BspTNI', 'AluBI', 'Van9116I', 'Psp0357II', 'BfaSII', 'Eco4174I', 'Lba2029III', 'CjeNIII', 'CjeI', 'PstNI', 'RruI', 'MunI', 'BssHII', 'Tru9I', 'AxyI', 'SgrBI', 'XagI', 'BcoDI', 'Bse1I', 'SwaI', 'ScoDS2II', 'NhaXI', 'CfrMH16VI', 'HdeZA17I', 'Hso63250IV', 'Eco81I', 'Bsu36I', 'Mly113I', 'Aor13HI', 'Eco72I', 'BglII', 'SacI', 'BstBAI', 'BspPI', 'CspCI', 'AleI', 'TpyTP2I', 'Bce10661III', 'Ecl35734I', 'Pru8113I', 'KpnNIH50I', 'CchIII', 'Sfr303I', 'PspCI', 'MseI', 'BsrGI', 'TatI', 'Ama87I', 'TseFI', 'PcsI', 'SnaBI', 'BccI', 'BfuI', 'Cfa8380I', 'Hca13221V', 'Sba460II', 'Nan12227I', 'Cin11811I', 'TspRI', 'PluTI', 'CspI', 'BstEII', 'Eco32I', 'BamHI', 'MauBI', 'AhlI', 'BssNAI', 'BslFI', 'BaeI', 'SurP32aII', 'Pme10899I', 'BbuB31II', 'DvuIII', 'KpnNH25III', 'BsbI', 'HpyF10VI', 'PmeI', 'MluI', 'BspTI', 'SspMI', 'SgeI', 'SdaI', 'Alw26I', 'AsuHPI', 'SfoI', 'Sau64037IV', 'MtuHN878II', 'Cdi13750III', 'GsuPI', 'SgrTI', 'BceSIV', 'CpoI', 'BsoBI', 'MaeI', 'AcsI', 'PflMI', 'DpnI', 'AscI', 'NsiI', 'BspFNI', 'BpiI', 'Sth20745III', 'PinP59III', 'Kpn9178I', 'Bbr7017III', 'DrdVIII', 'AquIV', 'SnaI', 'PdmI', 'MboI', 'BspDI', 'SlaI', 'SatI', 'RgaI', 'FriOI', 'EcoHI', 'PspOMII', 'RsaI', 'MspI7IV', 'Sau1803III', 'Cdi11397I', 'Gba708II', 'UbaF14I', 'EcoBLMcrX', 'Ksp22I', 'BlpI', 'AbsI', 'CviKI_1', 'ApaLI', 'Bst2UI', 'KpnI', 'Hpy188I', 'Tsp4CI', 'BbsI', 'Bsp68I', 'Bbr57III', 'DrdII', 'Pin17FIII', 'Kpn156V', 'Ssp6803IV', 'AquII', 'HgiEII', 'PceI', 'KasI', 'Bsp143I', 'SgrDI', 'Tth111I', 'PsyI', 'Psp124BI', 'Eam1105I', 'PsiI', 'Cco14983VI', 'FtnUV', 'Saf8902III', 'MspF392I', 'UbaF12I', 'GdiII', 'BinI', 'Tth111II', 'BstZ17I', 'AflII', 'Hsp92I', 'BanI', 'BstSCI', 'HaeII', 'EcoT38I', 'XmaIII', 'Psp03I', 'Acc36I', 'BshFI', 'AchA6III', 'Bau1417V', 'Dpi3084I', 'SpoDI', 'PfrJS15III', 'Kas9737III', 'AmaCSI', 'RleAI', 'NsbI', 'HpyCH4IV', 'Bsp119I', 'Sfr274I', 'TseI', 'PspEI', 'PaeI', 'DriI', 'PmlI', 'Mla10359I', 'Cco11437V', 'CjuI', 'Fnu11326IV', 'RspPBTS2III', 'UbaF9I', 'BbvII', 'TkoII', 'BstPAI', 'Acc65I', 'Hin6I', 'AvaI', 'BstENI', 'Cfr42I', 'DraIII', 'CfrI', 'McrI', 'TspGWI', 'BseJI', 'Spe19205IV', 'PfrJS12IV', 'AccX', 'Bag18758I', 'CspX1II', 'Jma19592II', 'EcoO157SI', 'MssI', 'HindIII', 'Bsp13I', 'RsaNI', 'StyD4I', 'Psp5II', 'Mph1103I', 'BstSLI', 'PspFI', 'NlaIV', 'Aco12261II', 'Cch467III', 'Rsp008V', 'Mch10819I', 'Fna13121I', 'RlaI', 'AceIII', 'TaqIII', 'FspBI', 'ApeKI', 'BstDEI', 'AsiSI', 'BstAPI', 'HgiAI', 'NmeAIII', 'BsaBI', 'Avi249I', 'CspBP25III', 'HspMHR1II', 'Sna507VIII', 'Pfl10783II', 'BscAI', 'BmeDI', 'MslI', 'HapII', 'BshVI', 'PsuI', 'SmlI', 'PflFI', 'Hsp92II', 'BstMCI', 'BseYI', 'SfaNI', 'MspA1I', 'Mba11I', 'Cbo67071IV', 'Fba202Z8II', 'RpaTI', 'Pfl1108I', 'SdeAI', 'XhoI', 'FauNDI', 'AflIII', 'BspT107I', 'AatII', 'BslI', 'Sth302II', 'BsiYI', 'SapI', 'MmeI', 'BmiI', 'SmaUMH5I', 'Pdu1735I', 'Asu14238IV', 'Csa9238II', 'HpyUM037X', 'Mox20I', 'EcoRI', 'BseX3I', 'PspLI', 'SexAI', 'PasI', 'FaeI', 'BsiHKAI', 'PciSI', 'Hpy8I', 'Rkr11038I', 'Lsp48III', 'Cba13II', 'EsaSSI', 'FinI', 'RpaBI', 'VspI', 'EclXI', 'BsiHKCI', 'BglI', 'Pfl8569I', 'VpaK11AI', 'BbrPI', 'PaqCI', 'HpyAV', 'PctI', 'Sma325I', 'AspSLV7III', 'Cpe13170II', 'HpyUM032XIII_mut1', 'Pcr308II', 'MlsI', 'EaeI', 'BseAI', 'PshBI', 'Sau96I', 'MteI', 'CfoI', 'Bsh1285I', 'HindII', 'LweI', 'BspACI', 'Lra68I', 'CalB3II', 'Eli8509II', 'Rgo13296IV', 'CjeP659IV', 'RlaII', 'CjePI', 'TaqI', 'Csp6I', 'BseDI', 'BaeGI', 'MstI', 'Sse232I', 'Sse8647I', 'Aor51HI', 'HgaI', 'Eco57I', 'LmnI', 'Sgr7807I', 'Pbu13063II', 'AspDUT2V', 'Cpe2837III', 'HpyLIM9XVI', 'KspAI', 'CviQI', 'Bsa29I', 'PinAI', 'PspPI', 'MaeIII', 'BstKTI', 'BseLI', 'LguI', 'AciI', 'HaeIII', 'Rer8036II', 'Lpn11417II', 'Bve1B23I', 'EcoNIH6II', 'BmgI', 'RceI', 'SpeI', 'CciI', 'Bse21I', 'AgsI', 'HaeI', 'SelI', 'SecI', 'AfaI', 'BtrI', 'Esp3I', 'BsrI', 'EciI', 'Seq11824I', 'Pal408I', 'Asp337I', 'Cla11845III', 'HpyG272XV', 'HpyCH4V', 'ClaI', 'BlnI', 'PciI', 'PpuMI', 'KflI', 'BstH2I', 'BlsI', 'TssI', 'FauI', 'BsaXI', 'EcoRV', 'Lmo911II', 'Rba2021I', 'Bsp3004IV', 'AbaPBA3II', 'EcoHSI', 'AbaUMB2I', 'SspD5I', 'Sau3AI', 'BstZI', 'Bpu1102I', 'BspLU11I', 'SanDI', 'MspGI', 'AjiI', 'EarI', 'BsgI', 'BtsIMutI', 'Acc16I', 'Asl11923II', 'CjeNV', 'HpyAXVIII', 'SenTFIV', 'Pae10662III', 'GlaI', 'Cfr9I', 'BfaI', 'PalAI', 'MvaI', 'Hpy188III', 'BspMAI', 'AlwNI', 'Eam1104I', 'ArsI', 'CviJI', 'BscGI', 'AbaB8342IV', 'Eco43896II', 'Pst14472I', 'Lme32I', 'Yru12986I', 'NmeA6CIII', 'ZraI', 'PspOMI', 'BstMBI', 'BmgT120I', 'Asi256I', 'EspI', 'BspMI', 'BseMII', 'BsuI', 'AhyYL17I', 'CjeFV', 'SenA1673III', 'Pac19842II', 'HpyAXVI_mut1', 'FaiI', 'BstYI', 'AvrII', 'NspV', 'Fsp4HI', 'FblI', 'AspLEI', 'AfiI', 'BveI', 'BsuRI', 'Lde4408II', 'BloAII', 'Aba6411II', 'Eco9035I', 'Pst145I', 'Xca85IV', 'MaqI', 'SseBI', 'NotI', 'BstAUI', 'XspI', 'Bme1390I', 'XceI', 'DraII', 'BsmBI', 'BpmI', 'BsmI', 'Sen6480IV', 'OspHL35III', 'Adh6U21I', 'Cje54107III', 'Hpy300XI', 'FnuDII', 'EgeI', 'BspHI', 'AsuII', 'NarI', 'EcoRII', 'EcoT14I', 'AccB7I', 'PabI', 'BstV1I', 'BstC8I', 'VpaSKIII', 'BkrAM31DI', 'BdaI', 'Eco9009II', 'PspMR102II', 'Lcr047II', 'DrdIV', 'SdeOSI', 'SmiMI', 'NgoMIV', 'BstACI', 'XapI', 'BisI', 'SstI', 'TaaI', 'AsuI', 'BsrBI', 'BsaI', 'AcuI', 'BseMI', 'Cgl13032II', 'Hpy99XIV_mut1', 'Sdy9603I', 'ObaBS10I', 'CviRI', 'EcoICRI', 'BseCI', 'Asp718I', 'MroNI', 'EcoNI', 'Eco130I', 'Van91I', 'Bst6I', 'BalI', 'Lbr124II', 'Bga514I', 'Eco4465II', 'VchE4II', 'PspAT13III', 'CstMI', 'NgoAVIII', 'NhoI', 'RseI', 'NcoI', 'BssMI', 'Vha464I', 'BciT130I', 'XmiI', 'SphI', 'SduI', 'BfuAI', 'GsaI', 'Bse3DI', 'XmnI', 'Sdy5370I', 'NpeUS61II', 'Cfupf3II', 'Hpy99XIII', 'MjaIV', 'Eco105I', 'BsaHI', 'AoxI', 'MreI', 'Cfr13I', 'Eco88I', 'SacII', 'TauI', 'BthCI', 'BspQI', 'BplI', 'AluI', 'TspARh3I', 'Pse18267I', 'Kro7512II', 'Bco11035III', 'Eco1836I', 'CdpI', 'PspN4I', 'MspI', 'BssAI', 'Tru1I', 'AsuC2I', 'VpaK11BI', 'SgfI', 'PkrI', 'SspI', 'BceAI', 'BmuI', 'Sbo46I', 'CfrMH13II', 'HdeNY26I', 'Nbr128II', 'HpyUM032XIII', 'BscXI', 'Eco47III', 'BclI', 'Alw44I', 'MflI', 'BstNI', 'Eco47I', 'PstI', 'XcmI', 'Bst1107I', 'Bso31I', 'BcgI', 'AccII', 'KpnNIH30III', 'Bce3081I', 'RdeGBIII', 'Ecl234I', 'Sxy1780I', 'PpiP13II', 'Ppu21I', 'MluCI', 'BspT104I', 'TasI', 'AccB1I', 'SmoI', 'SfaAI', 'MhlI', 'SmaI', 'BbvI', 'BciVI', 'SauMJ015III', 'Nal45188II', 'Cdu23823II', 'HbaII', 'AvaIII', 'Sth132I', 'Ecl136II', 'AseI', 'AcyI', 'MaeII', 'BssECI', 'CsiI', 'NspI', 'SfiI', 'BspLI', 'BseXI', 'AloI', 'SthSt3II', 'PliMI', 'Kpn9644II', 'BbuB31I', 'UcoMSI', 'DspS02II', 'BfiI', 'PmaCI', 'MfeI', 'BspEI', 'SspDI', 'SfcI', 'RigI', 'HpyCH4III', 'HgiCI', 'ScaI', 'SchI', 'AlwI', 'AbaSI', 'Sau5656II', 'MspSC27II', 'Cdi13746V', 'Gru56503II', 'UbaPI', 'Ksp632I', 'DinI', 'ApoI', 'AcoI', 'Kzo9I', 'BsaJI', 'BtgI', 'NlaIII', 'MwoI', 'BspANI', 'BmsI', 'Dde51507I', 'Kpn327I', 'SspJOR1II', 'NmeDI', 'Bbr7017II', 'DrdVI', 'PinP23II', 'AquIII', 'PdiI', 'Kpn2I', 'Bsp1407I', 'SgsI', 'Rsr2I', 'PvuI', 'Eco24I']

DNA = Seq

class FormattedSeq:
    lower: Incomplete
    data: Incomplete
    linear: Incomplete
    klass: Incomplete
    def __init__(self, seq, linear: bool = True) -> None: ...
    def __len__(self) -> int: ...
    def __eq__(self, other): ...
    def circularise(self) -> None: ...
    def linearise(self) -> None: ...
    def to_linear(self): ...
    def to_circular(self): ...
    def is_linear(self): ...
    def finditer(self, pattern, size): ...
    def __getitem__(self, i): ...

class RestrictionType(type):
    def __init__(cls, name: str = '', bases=(), dct: Incomplete | None = None) -> None: ...
    def __add__(cls, other): ...
    def __truediv__(cls, other): ...
    def __rtruediv__(cls, other): ...
    def __floordiv__(cls, other): ...
    def __rfloordiv__(cls, other): ...
    def __len__(cls) -> int: ...
    def __hash__(cls): ...
    def __eq__(cls, other): ...
    def __ne__(cls, other): ...
    def __rshift__(cls, other): ...
    def __mod__(cls, other): ...
    def __ge__(cls, other): ...
    def __gt__(cls, other): ...
    def __le__(cls, other): ...
    def __lt__(cls, other): ...

class AbstractCut(RestrictionType):
    @classmethod
    def search(cls, dna, linear: bool = True): ...
    @classmethod
    def all_suppliers(cls) -> None: ...
    @classmethod
    def is_equischizomer(cls, other): ...
    @classmethod
    def is_neoschizomer(cls, other): ...
    @classmethod
    def is_isoschizomer(cls, other): ...
    @classmethod
    def equischizomers(cls, batch: Incomplete | None = None): ...
    @classmethod
    def neoschizomers(cls, batch: Incomplete | None = None): ...
    @classmethod
    def isoschizomers(cls, batch: Incomplete | None = None): ...
    @classmethod
    def frequency(cls): ...

class NoCut(AbstractCut):
    @classmethod
    def cut_once(cls): ...
    @classmethod
    def cut_twice(cls): ...
    @classmethod
    def characteristic(cls): ...

class OneCut(AbstractCut):
    @classmethod
    def cut_once(cls): ...
    @classmethod
    def cut_twice(cls): ...
    @classmethod
    def characteristic(cls): ...

class TwoCuts(AbstractCut):
    @classmethod
    def cut_once(cls): ...
    @classmethod
    def cut_twice(cls): ...
    @classmethod
    def characteristic(cls): ...

class Meth_Dep(AbstractCut):
    @classmethod
    def is_methylable(cls): ...

class Meth_Undep(AbstractCut):
    @classmethod
    def is_methylable(cls): ...

class Palindromic(AbstractCut):
    @classmethod
    def is_palindromic(cls): ...

class NonPalindromic(AbstractCut):
    @classmethod
    def is_palindromic(cls): ...

class Unknown(AbstractCut):
    @classmethod
    def catalyse(cls, dna, linear: bool = True) -> None: ...
    catalyze = catalyse
    @classmethod
    def is_blunt(cls): ...
    @classmethod
    def is_5overhang(cls): ...
    @classmethod
    def is_3overhang(cls): ...
    @classmethod
    def overhang(cls): ...
    @classmethod
    def compatible_end(cls): ...

class Blunt(AbstractCut):
    @classmethod
    def catalyse(cls, dna, linear: bool = True): ...
    catalyze = catalyse
    @classmethod
    def is_blunt(cls): ...
    @classmethod
    def is_5overhang(cls): ...
    @classmethod
    def is_3overhang(cls): ...
    @classmethod
    def overhang(cls): ...
    @classmethod
    def compatible_end(cls, batch: Incomplete | None = None): ...

class Ov5(AbstractCut):
    @classmethod
    def catalyse(cls, dna, linear: bool = True): ...
    catalyze = catalyse
    @classmethod
    def is_blunt(cls): ...
    @classmethod
    def is_5overhang(cls): ...
    @classmethod
    def is_3overhang(cls): ...
    @classmethod
    def overhang(cls): ...
    @classmethod
    def compatible_end(cls, batch: Incomplete | None = None): ...

class Ov3(AbstractCut):
    @classmethod
    def catalyse(cls, dna, linear: bool = True): ...
    catalyze = catalyse
    @classmethod
    def is_blunt(cls): ...
    @classmethod
    def is_5overhang(cls): ...
    @classmethod
    def is_3overhang(cls): ...
    @classmethod
    def overhang(cls): ...
    @classmethod
    def compatible_end(cls, batch: Incomplete | None = None): ...

class Defined(AbstractCut):
    @classmethod
    def is_defined(cls): ...
    @classmethod
    def is_ambiguous(cls): ...
    @classmethod
    def is_unknown(cls): ...
    @classmethod
    def elucidate(cls): ...

class Ambiguous(AbstractCut):
    @classmethod
    def is_defined(cls): ...
    @classmethod
    def is_ambiguous(cls): ...
    @classmethod
    def is_unknown(cls): ...
    @classmethod
    def elucidate(cls): ...

class NotDefined(AbstractCut):
    @classmethod
    def is_defined(cls): ...
    @classmethod
    def is_ambiguous(cls): ...
    @classmethod
    def is_unknown(cls): ...
    @classmethod
    def elucidate(cls): ...

class Commercially_available(AbstractCut):
    @classmethod
    def suppliers(cls) -> None: ...
    @classmethod
    def supplier_list(cls): ...
    @classmethod
    def buffers(cls, supplier) -> None: ...
    @classmethod
    def is_comm(cls): ...

class Not_available(AbstractCut):
    @staticmethod
    def suppliers() -> None: ...
    @classmethod
    def supplier_list(cls): ...
    @classmethod
    def buffers(cls, supplier) -> None: ...
    @classmethod
    def is_comm(cls): ...

class RestrictionBatch(set):
    mapping: Incomplete
    already_mapped: Incomplete
    suppliers: Incomplete
    def __init__(self, first=(), suppliers=()) -> None: ...
    def __contains__(self, other) -> bool: ...
    def __div__(self, other): ...
    def __rdiv__(self, other): ...
    def __truediv__(self, other): ...
    def __rtruediv__(self, other): ...
    def get(self, enzyme, add: bool = False): ...
    def lambdasplit(self, func): ...
    def add_supplier(self, letter) -> None: ...
    def current_suppliers(self): ...
    def __iadd__(self, other): ...
    def __add__(self, other): ...
    def remove(self, other): ...
    def add(self, other): ...
    def add_nocheck(self, other): ...
    def format(self, y): ...
    def is_restriction(self, y): ...
    def split(self, *classes, **bool): ...
    def elements(self): ...
    def as_string(self): ...
    @classmethod
    def suppl_codes(cls): ...
    @classmethod
    def show_codes(cls) -> None: ...
    def search(self, dna, linear: bool = True): ...

class Analysis(RestrictionBatch, PrintFormat):
    rb: Incomplete
    sequence: Incomplete
    linear: Incomplete
    def __init__(self, restrictionbatch=..., sequence=..., linear: bool = True) -> None: ...
    def format_output(self, dct: Incomplete | None = None, title: str = '', s1: str = ''): ...
    def print_that(self, dct: Incomplete | None = None, title: str = '', s1: str = '') -> None: ...
    Cmodulo: Incomplete
    PrefWidth: Incomplete
    def change(self, **what) -> None: ...
    def full(self, linear: bool = True): ...
    def blunt(self, dct: Incomplete | None = None): ...
    def overhang5(self, dct: Incomplete | None = None): ...
    def overhang3(self, dct: Incomplete | None = None): ...
    def defined(self, dct: Incomplete | None = None): ...
    def with_sites(self, dct: Incomplete | None = None): ...
    def without_site(self, dct: Incomplete | None = None): ...
    def with_N_sites(self, N, dct: Incomplete | None = None): ...
    def with_number_list(self, list, dct: Incomplete | None = None): ...
    def with_name(self, names, dct: Incomplete | None = None): ...
    def with_site_size(self, site_size, dct: Incomplete | None = None): ...
    def only_between(self, start, end, dct: Incomplete | None = None): ...
    def between(self, start, end, dct: Incomplete | None = None): ...
    def show_only_between(self, start, end, dct: Incomplete | None = None): ...
    def only_outside(self, start, end, dct: Incomplete | None = None): ...
    def outside(self, start, end, dct: Incomplete | None = None): ...
    def do_not_cut(self, start, end, dct: Incomplete | None = None): ...

CommOnly: Incomplete
NonComm: Incomplete
AllEnzymes: Incomplete

# Names in __all__ with no definition:
#   AanI
#   AarI
#   AasI
#   AatII
#   Aba13301I
#   Aba6411II
#   AbaB8342IV
#   AbaCIII
#   AbaPBA3II
#   AbaSI
#   AbaUMB2I
#   Abr4036II
#   AbsI
#   Acc16I
#   Acc36I
#   Acc65I
#   Acc65V
#   AccB1I
#   AccB7I
#   AccBSI
#   AccI
#   AccII
#   AccIII
#   AccIX
#   AccX
#   AceIII
#   AchA6III
#   AciI
#   AclI
#   AclWI
#   Aco12261II
#   AcoI
#   AcoY31II
#   AcsI
#   AcuI
#   AcvI
#   AcyI
#   AdeI
#   Adh6U21I
#   AfaI
#   AfeI
#   AfiI
#   AflII
#   AflIII
#   AgeI
#   AgsI
#   AhaIII
#   AhdI
#   AhlI
#   AhyRBAHI
#   AhyYL17I
#   AjiI
#   AjnI
#   AjuI
#   AleI
#   AlfI
#   AloI
#   AluBI
#   AluI
#   Alw21I
#   Alw26I
#   Alw44I
#   AlwFI
#   AlwI
#   AlwNI
#   Ama87I
#   AmaCSI
#   Aod1I
#   Aor13HI
#   Aor51HI
#   AoxI
#   ApaBI
#   ApaI
#   ApaLI
#   ApeKI
#   ApoI
#   ApyPI
#   AquII
#   AquIII
#   AquIV
#   ArsI
#   AscI
#   AseI
#   Asi256I
#   AsiGI
#   AsiSI
#   Asl11923II
#   Asp103I
#   Asp114pII
#   Asp337I
#   Asp700I
#   Asp718I
#   AspA2I
#   AspAMDIV
#   AspBHI
#   AspDUT2V
#   AspJHL3II
#   AspLEI
#   AspNIH4III
#   AspS9I
#   AspSLV7III
#   Asu14238IV
#   AsuC2I
#   AsuHPI
#   AsuI
#   AsuII
#   AsuNHI
#   AteTI
#   AvaI
#   AvaII
#   AvaIII
#   Avi249I
#   AvrII
#   Awo1030IV
#   AxyI
#   BaeGI
#   BaeI
#   Bag18758I
#   BalI
#   BamHI
#   BanI
#   BanII
#   BanLI
#   BarI
#   Bau1417V
#   BauI
#   Bbr52II
#   Bbr57III
#   Bbr7017II
#   Bbr7017III
#   BbrPI
#   BbsI
#   BbuB31I
#   BbuB31II
#   Bbv12I
#   BbvCI
#   BbvI
#   BbvII
#   BccI
#   Bce10661III
#   Bce3081I
#   Bce83I
#   BceAI
#   BceSIV
#   BcefI
#   BcgI
#   BciT130I
#   BciVI
#   BclI
#   BcnI
#   Bco11035III
#   BcoDI
#   BcuI
#   BdaI
#   BetI
#   BfaI
#   BfaSII
#   BfiI
#   BfmI
#   BfoI
#   BfrI
#   BfuAI
#   BfuI
#   Bga514I
#   BglI
#   BglII
#   Bhe175II
#   BinI
#   BisI
#   BkrAM31DI
#   Ble402II
#   BlnI
#   BloAII
#   BlpI
#   BlsI
#   BmcAI
#   Bme1390I
#   Bme18I
#   BmeDI
#   BmeRI
#   BmeT110I
#   BmgBI
#   BmgI
#   BmgT120I
#   BmiI
#   BmrFI
#   BmrI
#   BmsI
#   BmtI
#   BmuI
#   BoxI
#   BpiI
#   BplI
#   BpmI
#   Bps6700III
#   Bpu10I
#   Bpu1102I
#   Bpu14I
#   BpuEI
#   BpuMI
#   Bsa29I
#   BsaAI
#   BsaBI
#   BsaHI
#   BsaI
#   BsaJI
#   BsaWI
#   BsaXI
#   BsbI
#   Bsc4I
#   BscAI
#   BscGI
#   BscXI
#   Bse118I
#   Bse1I
#   Bse21I
#   Bse3DI
#   Bse8I
#   BseAI
#   BseBI
#   BseCI
#   BseDI
#   BseGI
#   BseJI
#   BseLI
#   BseMI
#   BseMII
#   BseNI
#   BsePI
#   BseRI
#   BseSI
#   BseX3I
#   BseXI
#   BseYI
#   BsgI
#   Bsh1236I
#   Bsh1285I
#   BshFI
#   BshNI
#   BshTI
#   BshVI
#   BsiEI
#   BsiHKAI
#   BsiHKCI
#   BsiI
#   BsiSI
#   BsiWI
#   BsiYI
#   BslFI
#   BslI
#   BsmAI
#   BsmBI
#   BsmFI
#   BsmI
#   BsnI
#   Bso31I
#   BsoBI
#   Bsp119I
#   Bsp120I
#   Bsp1286I
#   Bsp13I
#   Bsp1407I
#   Bsp143I
#   Bsp1720I
#   Bsp19I
#   Bsp24I
#   Bsp3004IV
#   Bsp460III
#   Bsp68I
#   BspACI
#   BspANI
#   BspCNI
#   BspD6I
#   BspDI
#   BspEI
#   BspFNI
#   BspGI
#   BspHI
#   BspLI
#   BspLU11I
#   BspMAI
#   BspMI
#   BspMII
#   BspNCI
#   BspOI
#   BspPI
#   BspQI
#   BspT104I
#   BspT107I
#   BspTI
#   BspTNI
#   BsrBI
#   BsrDI
#   BsrFI
#   BsrGI
#   BsrI
#   BssAI
#   BssECI
#   BssHII
#   BssMI
#   BssNAI
#   BssNI
#   BssSI
#   BssT1I
#   Bst1107I
#   Bst2BI
#   Bst2UI
#   Bst4CI
#   Bst6I
#   BstACI
#   BstAFI
#   BstAPI
#   BstAUI
#   BstBAI
#   BstBI
#   BstC8I
#   BstDEI
#   BstDSI
#   BstEII
#   BstENI
#   BstF5I
#   BstFNI
#   BstH2I
#   BstHHI
#   BstKTI
#   BstMAI
#   BstMBI
#   BstMCI
#   BstMWI
#   BstNI
#   BstNSI
#   BstPAI
#   BstPI
#   BstSCI
#   BstSFI
#   BstSLI
#   BstSNI
#   BstUI
#   BstV1I
#   BstV2I
#   BstX2I
#   BstXI
#   BstYI
#   BstZ17I
#   BstZI
#   Bsu15I
#   Bsu36I
#   BsuI
#   BsuRI
#   BsuTUI
#   BtgI
#   BtgZI
#   BthCI
#   BtrI
#   BtsCI
#   BtsI
#   BtsIMutI
#   BtuMI
#   Bve1B23I
#   BveI
#   Cac8I
#   CaiI
#   Cal14237I
#   CalB3II
#   Cau10061II
#   CauII
#   Cba13II
#   Cba16038I
#   Cbo67071IV
#   CcaP7V
#   Cch467III
#   CchII
#   CchIII
#   CciI
#   CciNI
#   Cco11366VI
#   Cco11437V
#   Cco14983V
#   Cco14983VI
#   CcrNAIII
#   Cdi11397I
#   Cdi13746V
#   Cdi13750III
#   CdiI
#   CdpI
#   Cdu23823II
#   Cfa8380I
#   CfoI
#   Cfr10I
#   Cfr13I
#   Cfr42I
#   Cfr9I
#   CfrI
#   CfrMH13II
#   CfrMH16VI
#   Cfupf3II
#   Cgl13032I
#   Cgl13032II
#   ChaI
#   Cin11811I
#   Cje265V
#   Cje54107III
#   CjeFIII
#   CjeFV
#   CjeI
#   CjeNII
#   CjeNIII
#   CjeNV
#   CjeP659IV
#   CjePI
#   CjuI
#   CjuII
#   Cko11077IV
#   Cla11845III
#   ClaI
#   Cly7489II
#   Cpe10578V
#   Cpe13170II
#   Cpe2837III
#   CpoI
#   Cre7908I
#   Csa9238II
#   CseI
#   CsiI
#   Csp2014I
#   Csp6I
#   CspAI
#   CspBP25III
#   CspCI
#   CspI
#   CspL61I
#   CspX1II
#   CstMI
#   CviAII
#   CviJI
#   CviKI_1
#   CviQI
#   CviRI
#   Dde51507I
#   DdeI
#   DinI
#   Dpi3069I
#   Dpi3084I
#   Dpi3090II
#   DpnI
#   DpnII
#   DraI
#   DraII
#   DraIII
#   DraRI
#   DrdI
#   DrdII
#   DrdIV
#   DrdV
#   DrdVI
#   DrdVIII
#   DriI
#   DsaI
#   DseDI
#   DspS02II
#   DvuIII
#   EaeI
#   EagI
#   Eam1104I
#   Eam1105I
#   EarI
#   EciI
#   Ecl136II
#   Ecl234I
#   Ecl35734I
#   EclXI
#   Eco105I
#   Eco130I
#   Eco147I
#   Eco1836I
#   Eco24I
#   Eco31I
#   Eco32I
#   Eco4174I
#   Eco43896II
#   Eco4465II
#   Eco47I
#   Eco47III
#   Eco52I
#   Eco53kI
#   Eco57I
#   Eco57MI
#   Eco72I
#   Eco8164I
#   Eco81I
#   Eco88I
#   Eco9009II
#   Eco9020I
#   Eco9035I
#   Eco91I
#   Eco9699II
#   EcoBLMcrX
#   EcoE1140I
#   EcoHI
#   EcoHSI
#   EcoICRI
#   EcoMVII
#   EcoNI
#   EcoNIH6II
#   EcoO109I
#   EcoO157SI
#   EcoO65I
#   EcoRI
#   EcoRII
#   EcoRV
#   EcoT14I
#   EcoT22I
#   EcoT38I
#   EgeI
#   EheI
#   Ehi46392I
#   Eli8509II
#   ErhG4T10I
#   ErhI
#   EsaBC3I
#   EsaSSI
#   Esp3007I
#   Esp3I
#   EspI
#   FaeI
#   FaiI
#   FalI
#   FaqI
#   FatI
#   FauI
#   FauNDI
#   Fba202Z8II
#   FbaI
#   FblI
#   Fco1691IV
#   FinI
#   FmuI
#   Fna13121I
#   Fnu11326II
#   Fnu11326IV
#   Fnu4HI
#   FnuDII
#   FokI
#   FriOI
#   FseI
#   Fsp4HI
#   FspAI
#   FspBI
#   FspEI
#   FspI
#   FspPK15I
#   FtnUV
#   GauT27I
#   Gba708II
#   GdiII
#   GlaI
#   GluI
#   Gru56503II
#   GsaI
#   GsuI
#   GsuPI
#   HaeI
#   HaeII
#   HaeIII
#   HapII
#   HauII
#   HbaII
#   Hca13221V
#   HdeNY26I
#   HdeZA17I
#   HgaI
#   HgiAI
#   HgiCI
#   HgiEII
#   HgiJII
#   HhaI
#   Hin1I
#   Hin1II
#   Hin4I
#   Hin4II
#   Hin6I
#   HinP1I
#   HincII
#   HindII
#   HindIII
#   HinfI
#   HpaI
#   HpaII
#   HphI
#   Hpy166II
#   Hpy178III
#   Hpy188I
#   Hpy188III
#   Hpy300XI
#   Hpy8I
#   Hpy99I
#   Hpy99XIII
#   Hpy99XIV
#   Hpy99XIV_mut1
#   Hpy99XXII
#   HpyAS001VI
#   HpyAV
#   HpyAXIV
#   HpyAXVIII
#   HpyAXVI_mut1
#   HpyAXVI_mut2
#   HpyCH4III
#   HpyCH4IV
#   HpyCH4V
#   HpyF10VI
#   HpyF3I
#   HpyG272XV
#   HpyLIM6XII
#   HpyLIM9XVI
#   HpyPU007XIX
#   HpySE526I
#   HpyUM032XIII
#   HpyUM032XIII_mut1
#   HpyUM032XIV
#   HpyUM037X
#   Hso63250IV
#   Hso63373III
#   Hsp92I
#   Hsp92II
#   HspAI
#   HspMHR1II
#   Jma19592I
#   Jma19592II
#   Jsp2502II
#   Kas9737III
#   KasI
#   KflI
#   Kor51II
#   Kpn156V
#   Kpn2I
#   Kpn327I
#   Kpn9178I
#   Kpn9644II
#   KpnI
#   KpnNH25III
#   KpnNIH30III
#   KpnNIH50I
#   Kro7512II
#   KroI
#   KroNI
#   Ksp22I
#   Ksp632I
#   KspAI
#   KspI
#   Kzo9I
#   Lba2029III
#   Lbr124II
#   Lcr047I
#   Lcr047II
#   LcrJM4II
#   Lde4408II
#   LguI
#   LlaG50I
#   Lme32I
#   LmnI
#   Lmo370I
#   Lmo911II
#   Lpl1004II
#   Lpn11417II
#   Lpn12272I
#   LpnI
#   LpnPI
#   Lra68I
#   LsaDS4I
#   Lsp1109I
#   Lsp48III
#   Lsp6406VI
#   LweI
#   MabI
#   MaeI
#   MaeII
#   MaeIII
#   MalI
#   MaqI
#   MauBI
#   Mba11I
#   MbiI
#   MboI
#   MboII
#   McaTI
#   Mch10819I
#   Mch946II
#   McrI
#   MfeI
#   MflI
#   MhlI
#   MjaIV
#   MkaDII
#   Mla10359I
#   MlsI
#   Mlu211III
#   MluCI
#   MluI
#   MluNI
#   Mly113I
#   MlyI
#   MmeI
#   MnlI
#   Mox20I
#   Mph1103I
#   MreI
#   MroI
#   MroNI
#   MroXI
#   MscI
#   MseI
#   MslI
#   Msp20I
#   MspA1I
#   MspCI
#   MspF392I
#   MspGI
#   MspI
#   MspI7II
#   MspI7IV
#   MspJI
#   MspR9I
#   MspSC27II
#   MssI
#   MstI
#   MteI
#   MtuHN878II
#   MunI
#   Mva1269I
#   MvaI
#   MvnI
#   MwoI
#   NaeI
#   Nal45188II
#   Nan12227I
#   NarI
#   Nbr128II
#   NciI
#   NcoI
#   NdeI
#   NdeII
#   NgoAVII
#   NgoAVIII
#   NgoMIV
#   NhaXI
#   NheI
#   NhoI
#   NlaCI
#   NlaIII
#   NlaIV
#   Nli3877I
#   NmeA6CIII
#   NmeAIII
#   NmeDI
#   NmuCI
#   NotI
#   NpeUS61II
#   NruI
#   NsbI
#   NsiI
#   NspBII
#   NspES21II
#   NspI
#   NspV
#   ObaBS10I
#   OgrI
#   OliI
#   OspHL35III
#   PabI
#   Pac19842II
#   PacI
#   PacIII
#   Pae10662III
#   Pae8506I
#   PaeI
#   PaePA99III
#   PaeR7I
#   PagI
#   Pal408I
#   PalAI
#   PaqCI
#   PasI
#   PauI
#   Pba2294I
#   Pbu13063II
#   PcaII
#   PceI
#   PciI
#   PciSI
#   Pcr308II
#   PcsI
#   PctI
#   Pdi8503III
#   PdiI
#   PdmI
#   Pdu1735I
#   PenI
#   PfeI
#   Pfl10783II
#   Pfl1108I
#   Pfl23II
#   Pfl3756II
#   Pfl8569I
#   PflFI
#   PflMI
#   PflPt14I
#   PfoI
#   PfrJS12IV
#   PfrJS12V
#   PfrJS15III
#   PgaP73III
#   Pin17FIII
#   PinAI
#   PinP23II
#   PinP59III
#   PkrI
#   PlaDI
#   Ple19I
#   PleI
#   PliMI
#   PluTI
#   PmaCI
#   Pme10899I
#   PmeI
#   PmlI
#   PpiI
#   PpiP13II
#   PpsI
#   Ppu10I
#   Ppu21I
#   PpuMI
#   Pru8113I
#   PscI
#   Pse18267I
#   PshAI
#   PshBI
#   PsiI
#   Psp0357II
#   Psp03I
#   Psp124BI
#   Psp1406I
#   Psp5II
#   Psp6I
#   PspAT13III
#   PspCI
#   PspD7DII
#   PspEI
#   PspFI
#   PspGI
#   PspLI
#   PspMR102II
#   PspN4I
#   PspOMI
#   PspOMII
#   PspPI
#   PspPPI
#   PspPRI
#   PspR84I
#   PspXI
#   PsrI
#   PssI
#   Pst14472I
#   Pst145I
#   Pst273I
#   PstI
#   PstNI
#   PsuGI
#   PsuI
#   PsyI
#   PteI
#   PvuI
#   PvuII
#   Ran11014IV
#   Rba2021I
#   RceI
#   RdeGBI
#   RdeGBII
#   RdeGBIII
#   Rer8036II
#   RflFIII
#   RgaI
#   Rgo13296IV
#   Rho5650I
#   RigI
#   Rkr11038I
#   RlaI
#   RlaII
#   RleAI
#   Rmu369III
#   RpaB5I
#   RpaBI
#   RpaI
#   RpaTI
#   RruI
#   RsaI
#   RsaNI
#   RseI
#   Rsp008IV
#   Rsp008V
#   Rsp531II
#   RspPBTS2III
#   Rsr2I
#   RsrII
#   Rtr1953I
#   SacI
#   SacII
#   Saf8902III
#   Sag901I
#   SalI
#   SanDI
#   SapI
#   SaqAI
#   SatI
#   Sau1803III
#   Sau3AI
#   Sau5656II
#   Sau64037IV
#   Sau96I
#   SauI
#   SauMJ015III
#   Sba460II
#   SbfI
#   Sbo46I
#   ScaI
#   SchI
#   SciI
#   ScoDS2II
#   ScrFI
#   SdaI
#   SdeAI
#   SdeOSI
#   SduI
#   Sdy5370I
#   Sdy7136I
#   Sdy9603I
#   SecI
#   SelI
#   Sen17963III
#   Sen5794III
#   Sen6480IV
#   SenA1673III
#   SenSARA26III
#   SenTFIV
#   Sep11964I
#   Seq11824I
#   SetI
#   SexAI
#   SfaAI
#   SfaNI
#   SfcI
#   SfeI
#   SfiI
#   Sfl13829III
#   SfoI
#   Sfr274I
#   Sfr303I
#   SfuI
#   SgeI
#   SgfI
#   Sgr7807I
#   SgrAI
#   SgrAII
#   SgrBI
#   SgrDI
#   SgrTI
#   SgsI
#   SimI
#   SinI
#   SlaI
#   Sma10259II
#   Sma325I
#   SmaI
#   SmaUMH5I
#   SmaUMH8I
#   SmiI
#   SmiMI
#   SmlI
#   SmoI
#   Sna507VIII
#   SnaBI
#   SnaI
#   Sno506I
#   Spe19205IV
#   SpeI
#   SphI
#   SplI
#   SpnRII
#   SpoDI
#   SrfI
#   Sse232I
#   Sse8387I
#   Sse8647I
#   Sse9I
#   SseBI
#   SsiI
#   Ssp6803IV
#   Ssp714II
#   SspD5I
#   SspDI
#   SspI
#   SspJOR1II
#   SspMI
#   SstE37I
#   SstI
#   Sth132I
#   Sth20745III
#   Sth302II
#   SthSt3II
#   StsI
#   StuI
#   StyD4I
#   StyI
#   SurP32aII
#   SwaI
#   Sxy1780I
#   TaaI
#   TagI
#   TaiI
#   TaqI
#   TaqII
#   TaqIII
#   TasI
#   TatI
#   TauI
#   TfiI
#   TkoI
#   TkoII
#   TpyTP2I
#   Tru1I
#   Tru9I
#   TscAI
#   TseFI
#   TseI
#   TsoI
#   Tsp45I
#   Tsp4CI
#   TspARh3I
#   TspDTI
#   TspEI
#   TspGWI
#   TspMI
#   TspRI
#   TssI
#   TstI
#   TsuI
#   Tth111I
#   Tth111II
#   UbaF11I
#   UbaF12I
#   UbaF13I
#   UbaF14I
#   UbaF9I
#   UbaPI
#   UcoMSI
#   UnbI
#   UpaP162I
#   Van9116I
#   Van91I
#   VchE4II
#   Vdi96II
#   Vha464I
#   VneI
#   VpaK11AI
#   VpaK11BI
#   VpaSKIII
#   VspI
#   Vtu19109I
#   WviI
#   XagI
#   XapI
#   XbaI
#   Xca85IV
#   XceI
#   XcmI
#   XhoI
#   XhoII
#   XmaI
#   XmaIII
#   XmaJI
#   XmiI
#   XmnI
#   XspI
#   YkrI
#   Yps3606I
#   Yru12986I
#   ZraI
#   ZrmI
#   Zsp2I
