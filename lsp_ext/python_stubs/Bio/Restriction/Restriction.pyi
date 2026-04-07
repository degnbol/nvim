from Bio.Restriction.PrintFormat import PrintFormat
from Bio.Seq import Seq
from _typeshed import Incomplete

__all__ = ['FormattedSeq', 'Analysis', 'RestrictionBatch', 'AllEnzymes', 'CommOnly', 'NonComm', 'MspF392I', 'BspLU11I', 'Hpy99XIV_mut1', 'MslI', 'HapII', 'FalI', 'BauI', 'AlwI', 'BpiI', 'SnaBI', 'FbaI', 'XmaI', 'ApaI', 'FinI', 'Cpe10578V', 'Sau5656II', 'MhlI', 'Bbr57III', 'PacIII', 'HpyUM037X', 'AcoY31II', 'YkrI', 'VpaK11AI', 'DspS02II', 'CjuI', 'SenTFIV', 'DraRI', 'HgiAI', 'Bps6700III', 'SsiI', 'PflPt14I', 'LweI', 'AlfI', 'SmiMI', 'BspCNI', 'Aor51HI', 'NcoI', 'XmiI', 'AsiGI', 'BtsCI', 'BsaWI', 'Eco31I', 'FauI', 'BciT130I', 'KpnNIH50I', 'SimI', 'Aba13301I', 'EcoHSI', 'BspEI', 'SpoDI', 'RpaB5I', 'Cco11366VI', 'PspAT13III', 'StyI', 'SgrTI', 'HpyUM032XIII', 'Sth302II', 'Lpn12272I', 'BspMI', 'FtnUV', 'BstZI', 'NaeI', 'Vtu19109I', 'BmuI', 'Cgl13032II', 'BspMII', 'BbsI', 'BspLI', 'BseX3I', 'BstMCI', 'Mox20I', 'PspLI', 'AccI', 'XbaI', 'SfoI', 'EcoRI', 'Rkr11038I', 'AspNIH4III', 'MspSC27II', 'SplI', 'AccX', 'HpyAXIV', 'PsuGI', 'CauII', 'Csa9238II', 'Sba460II', 'ApyPI', 'BbuB31I', 'Pae10662III', 'BstMAI', 'PceI', 'Jma19592I', 'BbrPI', 'SduI', 'MspI', 'VpaK11BI', 'BglII', 'ArsI', 'Eam1104I', 'Aor13HI', 'BtsI', 'BtgZI', 'PspN4I', 'Ecl35734I', 'NhoI', 'MaqI', 'Nli3877I', 'Sfl13829III', 'Bsp3004IV', 'PfrJS15III', 'Lbr124II', 'PteI', 'AvaIII', 'AbaPBA3II', 'BarI', 'Ehi46392I', 'SspJOR1II', 'TaqIII', 'Cco14983VI', 'PspR84I', 'BscXI', 'AccIII', 'Sau96I', 'Acc36I', 'BsiHKAI', 'RsaI', 'PshBI', 'CspAI', 'MlsI', 'Eco72I', 'BspANI', 'BseAI', 'Vha464I', 'CfoI', 'Gru56503II', 'Yru12986I', 'CjeFIII', 'XmaIII', 'Rsp008IV', 'Bme18I', 'CaiI', 'MbiI', 'BseCI', 'Asu14238IV', 'Nan12227I', 'HpyAXVIII', 'AspJHL3II', 'EsaBC3I', 'TspDTI', 'EspI', 'CspL61I', 'Sdy5370I', 'MluCI', 'BamHI', 'Ppu21I', 'TaaI', 'BpuEI', 'PluTI', 'CspI', 'AhlI', 'BsrDI', 'BssT1I', 'Bce10661III', 'BveI', 'BstEII', 'Pba2294I', 'Kas9737III', 'PalAI', 'BoxI', 'BcgI', 'Eco4465II', 'BsiSI', 'UbaPI', 'Sma325I', 'PlaDI', 'Tsp4CI', 'Cal14237I', 'MseI', 'PinP23II', 'PspGI', 'AceIII', 'RdeGBIII', 'Cre7908I', 'Ksp632I', 'LcrJM4II', 'BceAI', 'BstFNI', 'McaTI', 'TspGWI', 'CciNI', 'Bsa29I', 'BssMI', 'BstKTI', 'BseBI', 'KspAI', 'BsnI', 'PsiI', 'FspI', 'MteI', 'SurP32aII', 'Cdi13746V', 'Sse9I', 'TsoI', 'Pst14472I', 'BlsI', 'AhyRBAHI', 'AflII', 'Pfl3756II', 'Mch946II', 'XmaJI', 'HgaI', 'BmgBI', 'EcoICRI', 'Hca13221V', 'EclXI', 'AlwFI', 'BsmI', 'CjeNV', 'AcoI', 'RspPBTS2III', 'MboII', 'Awo1030IV', 'BsrFI', 'NpeUS61II', 'AcsI', 'NsiI', 'BseNI', 'MfeI', 'AscI', 'CpoI', 'PmaCI', 'BmrI', 'DpnI', 'BstV1I', 'AjnI', 'MauBI', 'SgfI', 'HpyLIM6XII', 'BmeDI', 'SauI', 'Dpi3084I', 'BfaI', 'Sen5794III', 'BstENI', 'CchII', 'Van91I', 'Bga514I', 'HinP1I', 'UbaF13I', 'Pcr308II', 'Kpn327I', 'Pfl23II', 'UcoMSI', 'CjeI', 'Eco9020I', 'AluI', 'SmaUMH8I', 'RceI', 'ClaI', 'NmeAIII', 'BpuMI', 'HpyCH4V', 'PpuMI', 'Bsh1236I', 'PmlI', 'PciI', 'BsuTUI', 'BlnI', 'Cba13II', 'Pme10899I', 'Lsp48III', 'BcefI', 'LpnI', 'Lme32I', 'Hin4I', 'Fco1691IV', 'AccBSI', 'TspARh3I', 'Cfa8380I', 'RdeGBI', 'BpmI', 'Asl11923II', 'ApoI', 'Mla10359I', 'PleI', 'Asi256I', 'EheI', 'Hpy99XIII', 'SgsI', 'Bpu10I', 'PvuI', 'BstUI', 'HpyCH4III', 'Bsp1407I', 'AbsI', 'Rsr2I', 'BseGI', 'Kpn2I', 'Cly7489II', 'CjeP659IV', 'Sag901I', 'Bau1417V', 'OgrI', 'UbaF11I', 'HpyUM032XIII_mut1', 'AcvI', 'Hin4II', 'Sse8647I', 'DrdVI', 'Bpu14I', 'SenA1673III', 'CjeNIII', 'BsiYI', 'Ble402II', 'HpyCH4IV', 'Bst2BI', 'WviI', 'LguI', 'SalI', 'Hpy188III', 'Cfr9I', 'MmeI', 'ZrmI', 'MvaI', 'PpsI', 'GlaI', 'Bsc4I', 'SapI', 'RseI', 'KpnNH25III', 'BsiI', 'HauII', 'SdeOSI', 'Eco43896II', 'Spe19205IV', 'RpaI', 'CcaP7V', 'NgoAVII', 'Pse18267I', 'EcoBLMcrX', 'Sma10259II', 'Pfl8569I', 'Lpl1004II', 'BsmBI', 'DinI', 'Fnu11326IV', 'MscI', 'Vdi96II', 'BciVI', 'Cfupf3II', 'BstZ17I', 'BspQI', 'BcoDI', 'BsaAI', 'Hsp92I', 'BanI', 'Hpy188I', 'BseRI', 'Bsp120I', 'Ple19I', 'AsuC2I', 'Rgo13296IV', 'AspAMDIV', 'MspI7II', 'BscAI', 'Ppu10I', 'Hpy99XXII', 'PenI', 'BspACI', 'EciI', 'TsuI', 'Cpe13170II', 'AoxI', 'Sau64037IV', 'Bbr7017II', 'CdiI', 'Pac19842II', 'BspTNI', 'NsbI', 'Hso63373III', 'AspLEI', 'PctI', 'PaeR7I', 'FblI', 'BstYI', 'StuI', 'HpyAV', 'MspJI', 'FaiI', 'Fsp4HI', 'DvuIII', 'Sep11964I', 'Dde51507I', 'DrdIV', 'BscGI', 'HgiJII', 'PfrJS12IV', 'BdaI', 'MlyI', 'Kro7512II', 'AbaB8342IV', 'AjuI', 'EcoMVII', 'BsuRI', 'Ssp714II', 'SdeAI', 'Cco11437V', 'PspOMI', 'PspD7DII', 'Sth132I', 'MvnI', 'AvaI', 'BstPAI', 'Bso31I', 'XmnI', 'BccI', 'Hin6I', 'AluBI', 'EcoT38I', 'Cfr42I', 'Tru1I', 'Eco32I', 'Lra68I', 'GauT27I', 'NlaIV', 'Xca85IV', 'Bse1I', 'PspPPI', 'Cje265V', 'CfrI', 'Rmu369III', 'AciI', 'AspSLV7III', 'BsaHI', 'PspOMII', 'MtuHN878II', 'Sse232I', 'Acc65V', 'HpyAXVI_mut1', 'RlaI', 'Mva1269I', 'DraII', 'Csp2014I', 'BstAFI', 'LmnI', 'NheI', 'BspHI', 'SrfI', 'Eco57I', 'EgeI', 'LpnPI', 'Bsp1720I', 'TaiI', 'EaeI', 'Sbo46I', 'BbuB31II', 'PaePA99III', 'AquII', 'BsoBI', 'PdiI', 'ChaI', 'Jma19592II', 'NspV', 'BmcAI', 'Eco1836I', 'BshTI', 'SnaI', 'Sgr7807I', 'NlaCI', 'Psp03I', 'BspNCI', 'MluI', 'PgaP73III', 'Asp114pII', 'ZraI', 'Lcr047I', 'BbvI', 'Msp20I', 'BstBAI', 'BseXI', 'SspI', 'Alw26I', 'SbfI', 'AleI', 'BsaXI', 'AsiSI', 'AloI', 'Bst1107I', 'BssAI', 'DraI', 'MaeIII', 'Eli8509II', 'Sth20745III', 'CcrNAIII', 'SgrAI', 'TkoI', 'Pst145I', 'Tsp45I', 'AlwNI', 'Acc65I', 'Lsp6406VI', 'GsuPI', 'Csp6I', 'PshAI', 'AasI', 'CjeFV', 'Rsp008V', 'HphI', 'AteTI', 'BsiWI', 'MroNI', 'Asp718I', 'AccB7I', 'Sse8387I', 'BfmI', 'VneI', 'FspEI', 'Esp3I', 'BsrI', 'HindII', 'Nbr128II', 'SfaAI', 'Asp337I', 'HpyAS001VI', 'MaeI', 'FnuDII', 'Hpy178III', 'CspX1II', 'AvrII', 'Sdy7136I', 'BstDEI', 'BfiI', 'TauI', 'Bco11035III', 'Hso63250IV', 'Pbu13063II', 'CseI', 'Kor51II', 'PauI', 'CspCI', 'Eco8164I', 'PfeI', 'SinI', 'PspXI', 'BstMWI', 'Hin1II', 'FatI', 'BplI', 'BmsI', 'MroXI', 'BssNAI', 'AccII', 'Eco91I', 'CalB3II', 'PspPI', 'PinP59III', 'TspRI', 'EsaSSI', 'AquIV', 'AspBHI', 'AhaIII', 'Lde4408II', 'SfuI', 'Bsp24I', 'PabI', 'Esp3007I', 'BssNI', 'Sxy1780I', 'Tth111II', 'Cdi13750III', 'Ran11014IV', 'AcuI', 'EcoHI', 'AhyYL17I', 'Mch10819I', 'Lsp1109I', 'Eco53kI', 'MreI', 'SacII', 'Eco88I', 'Cfr13I', 'Eco105I', 'Eco52I', 'AfeI', 'HaeIII', 'EarI', 'FaqI', 'Hsp92II', 'SmoI', 'HdeNY26I', 'Cko11077IV', 'BmgI', 'Rtr1953I', 'BmrFI', 'Eco24I', 'Bag18758I', 'NspES21II', 'ApeKI', 'Cin11811I', 'SgrBI', 'HpyLIM9XVI', 'Eco57MI', 'SecI', 'AsuI', 'Dpi3090II', 'BfrI', 'Sen6480IV', 'CchIII', 'Bhe175II', 'HindIII', 'BsiHKCI', 'EcoT22I', 'NmuCI', 'VspI', 'EagI', 'AclWI', 'MluNI', 'BspFNI', 'ScaI', 'BsePI', 'DdeI', 'BanII', 'Bsp1286I', 'PinAI', 'Kpn9178I', 'KpnI', 'NgoAVIII', 'Eco9035I', 'BstXI', 'Sna507VIII', 'EcoT14I', 'RdeGBII', 'Cba16038I', 'SchI', 'NdeI', 'PpiP13II', 'ScrFI', 'BinI', 'BthCI', 'Pdi8503III', 'MstI', 'Lmo370I', 'BsaI', 'BtuMI', 'Fna13121I', 'HpaI', 'PstNI', 'AbaSI', 'BsrGI', 'TatI', 'PstI', 'Eco47I', 'PspCI', 'BstNI', 'MflI', 'Eco47III', 'PflFI', 'AccB1I', 'CfrMH13II', 'Bsh1285I', 'Rer8036II', 'Asp103I', 'Mlu211III', 'BetI', 'PacI', 'Hpy99XIV', 'FspBI', 'SmaI', 'CjuII', 'BtsIMutI', 'XagI', 'Cpe2837III', 'Alw44I', 'Sau1803III', 'Bbr52II', 'OspHL35III', 'BspPI', 'MssI', 'BseDI', 'BaeGI', 'TaqI', 'DpnII', 'MalI', 'AarI', 'Bsp68I', 'PvuII', 'Bse118I', 'AfaI', 'MroI', 'DrdVIII', 'RleAI', 'UnbI', 'SenSARA26III', 'Bst2UI', 'CstMI', 'BloAII', 'FmuI', 'HspAI', 'PspFI', 'Pfl10783II', 'EcoRII', 'DraIII', 'BtrI', 'SmiI', 'BslI', 'KpnNIH30III', 'Psp1406I', 'GdiII', 'Aba6411II', 'TstI', 'EcoE1140I', 'BstC8I', 'SpnRII', 'RpaBI', 'Cch467III', 'NotI', 'Psp0357II', 'SspMI', 'PcsI', 'BspTI', 'SdaI', 'SgeI', 'NspI', 'PmeI', 'Ecl136II', 'AanI', 'CviJI', 'SmlI', 'BsmFI', 'TasI', 'Lpn11417II', 'SciI', 'BstH2I', 'FspPK15I', 'BstX2I', 'MspA1I', 'VpaSKIII', 'BfuI', 'Psp6I', 'Cgl13032I', 'UbaF14I', 'Rho5650I', 'BsgI', 'AxyI', 'Bst4CI', 'AspDUT2V', 'BclI', 'MspI7IV', 'StsI', 'SelI', 'Hpy166II', 'Hpy300XI', 'AccIX', 'Pfl1108I', 'GsuI', 'BseSI', 'MabI', 'BstHHI', 'AgsI', 'Bse21I', 'SpeI', 'CviAII', 'KroNI', 'BshFI', 'TaqII', 'BshNI', 'PkrI', 'SauMJ015III', 'CviQI', 'Bbr7017III', 'AmaCSI', 'Pae8506I', 'SspD5I', 'Bst6I', 'BlpI', 'OliI', 'Zsp2I', 'HspMHR1II', 'NarI', 'AjiI', 'Asp700I', 'AatII', 'BceSIV', 'Ecl234I', 'Seq11824I', 'HgiEII', 'CsiI', 'DrdV', 'Bsp460III', 'McrI', 'MboI', 'PfrJS12V', 'PciSI', 'Hpy99I', 'SseBI', 'Lba2029III', 'BtgI', 'BsaJI', 'Kzo9I', 'SfiI', 'RgaI', 'SatI', 'SlaI', 'HpyF10VI', 'BspDI', 'BseMI', 'PsuI', 'BspT104I', 'SacI', 'AbaCIII', 'EcoNIH6II', 'Cac8I', 'Ssp6803IV', 'KflI', 'SstE37I', 'Cco14983V', 'PspMR102II', 'TfiI', 'AfiI', 'UbaF12I', 'LsaDS4I', 'Tru9I', 'BstNSI', 'Gba708II', 'CciI', 'Yps3606I', 'Cje54107III', 'XhoII', 'RpaTI', 'BstF5I', 'BspOI', 'HpyF3I', 'NciI', 'BmeRI', 'BseJI', 'NruI', 'SfaNI', 'MnlI', 'Sau3AI', 'Cfr10I', 'BstSLI', 'Nal45188II', 'TspEI', 'RigI', 'AchA6III', 'HpyAXVI_mut2', 'Ksp22I', 'CviRI', 'DsaI', 'CspBP25III', 'AsuII', 'ScoDS2II', 'BspT107I', 'AquIII', 'SetI', 'Bce3081I', 'Pal408I', 'Van9116I', 'BstV2I', 'BssECI', 'PdmI', 'MspGI', 'Jsp2502II', 'BmiI', 'BaeI', 'Eco4174I', 'HhaI', 'BstSFI', 'BcnI', 'KroI', 'MwoI', 'Psp124BI', 'PsyI', 'Bse3DI', 'KasI', 'AgeI', 'BshVI', 'Eco81I', 'SgrAII', 'NmeA6CIII', 'Bve1B23I', 'Pin17FIII', 'PssI', 'PflMI', 'Lcr047II', 'SaqAI', 'BssSI', 'PsrI', 'UbaF9I', 'ErhG4T10I', 'BssHII', 'EcoRV', 'SthSt3II', 'MspR9I', 'TkoII', 'Cdi11397I', 'Pst273I', 'GsaI', 'Tth111I', 'Adh6U21I', 'AclI', 'Mba11I', 'FokI', 'BstMBI', 'BmgT120I', 'BfoI', 'Eco147I', 'GluI', 'HinfI', 'PagI', 'Bbv12I', 'FspAI', 'BsaBI', 'FaeI', 'HbaII', 'AbaUMB2I', 'SfcI', 'CjeNII', 'Rsp531II', 'BmeT110I', 'DseDI', 'Avi249I', 'NhaXI', 'AflIII', 'Sfr303I', 'HpyG272XV', 'MaeII', 'UpaP162I', 'Bce83I', 'SanDI', 'Dpi3069I', 'BcuI', 'TssI', 'Sdy9603I', 'BstDSI', 'BsbI', 'BfaSII', 'Bsp119I', 'PspEI', 'TseI', 'PaeI', 'Sfr274I', 'Eam1105I', 'FseI', 'BstPI', 'AvaII', 'HpySE526I', 'Bsu36I', 'PcaII', 'AhdI', 'Kpn156V', 'Bse8I', 'HaeII', 'Eco9009II', 'Bsp19I', 'TagI', 'SmaUMH5I', 'BbvCI', 'Eco130I', 'Abr4036II', 'PspPRI', 'Cau10061II', 'MunI', 'PliMI', 'RsrII', 'BbvII', 'XcmI', 'HpyUM032XIV', 'HaeI', 'LlaG50I', 'BfuAI', 'SgrDI', 'BstSNI', 'CjePI', 'Fba202Z8II', 'XceI', 'Alw21I', 'NdeII', 'Fnu4HI', 'HincII', 'AsuNHI', 'Bme1390I', 'XspI', 'BstAUI', 'TspMI', 'PasI', 'Ama87I', 'TpyTP2I', 'Cdu23823II', 'Rba2021I', 'HgiCI', 'NmeDI', 'BseLI', 'Aod1I', 'ApaLI', 'MkaDII', 'PaqCI', 'KspI', 'BsrBI', 'HdeZA17I', 'FauNDI', 'BspGI', 'BsuI', 'TseFI', 'Cla11845III', 'AcyI', 'Saf8902III', 'Bpu1102I', 'FriOI', 'BanLI', 'Bsu15I', 'ObaBS10I', 'AspS9I', 'BslFI', 'Mph1103I', 'Psp5II', 'StyD4I', 'Hin1I', 'DrdI', 'RsaNI', 'DriI', 'Bsp13I', 'SwaI', 'BmtI', 'Mly113I', 'Acc16I', 'SphI', 'HpyPU007XIX', 'EcoO157SI', 'DrdII', 'SfeI', 'Sen17963III', 'BstSCI', 'CdpI', 'BkrAM31DI', 'ApaBI', 'HpaII', 'BseYI', 'Pdu1735I', 'EcoNI', 'BstAPI', 'RruI', 'BglI', 'Kpn9644II', 'PscI', 'NlaIII', 'PpiI', 'Eco9699II', 'Bsp143I', 'BalI', 'Sno506I', 'ErhI', 'RlaII', 'Cbo67071IV', 'NgoMIV', 'EcoO65I', 'EcoO109I', 'MspCI', 'AdeI', 'AspA2I', 'SstI', 'BisI', 'XapI', 'TscAI', 'BstACI', 'SexAI', 'Pru8113I', 'BspD6I', 'Lmo911II', 'NspBII', 'BsmAI', 'SspDI', 'CviKI_1', 'Aco12261II', 'BspMAI', 'Fnu11326II', 'BstBI', 'Hpy8I', 'VchE4II', 'AsuHPI', 'PfoI', 'CfrMH16VI', 'MjaIV', 'XhoI', 'RflFIII', 'BseMII', 'BsiEI', 'AseI']

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
    def __init__(cls, name: str = '', bases=(), dct=None) -> None: ...
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
    def equischizomers(cls, batch=None): ...
    @classmethod
    def neoschizomers(cls, batch=None): ...
    @classmethod
    def isoschizomers(cls, batch=None): ...
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
    def compatible_end(cls, batch=None): ...

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
    def compatible_end(cls, batch=None): ...

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
    def compatible_end(cls, batch=None): ...

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
    def format_output(self, dct=None, title: str = '', s1: str = ''): ...
    def print_that(self, dct=None, title: str = '', s1: str = '') -> None: ...
    Cmodulo: Incomplete
    PrefWidth: Incomplete
    def change(self, **what) -> None: ...
    def full(self, linear: bool = True): ...
    def blunt(self, dct=None): ...
    def overhang5(self, dct=None): ...
    def overhang3(self, dct=None): ...
    def defined(self, dct=None): ...
    def with_sites(self, dct=None): ...
    def without_site(self, dct=None): ...
    def with_N_sites(self, N, dct=None): ...
    def with_number_list(self, list, dct=None): ...
    def with_name(self, names, dct=None): ...
    def with_site_size(self, site_size, dct=None): ...
    def only_between(self, start, end, dct=None): ...
    def between(self, start, end, dct=None): ...
    def show_only_between(self, start, end, dct=None): ...
    def only_outside(self, start, end, dct=None): ...
    def outside(self, start, end, dct=None): ...
    def do_not_cut(self, start, end, dct=None): ...

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
