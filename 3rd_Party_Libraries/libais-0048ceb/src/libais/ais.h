// -*- c++ -*-

// TODO(schwehr): Create an archive of messages that do not decode.

#ifndef LIBAIS_AIS_H_
#define LIBAIS_AIS_H_

#include <array>
#include <bitset>
#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>

using std::bitset;
using std::ostream;
using std::string;
using std::vector;

#define LIBAIS_VERSION_MAJOR 0
#define LIBAIS_VERSION_MINOR 15

extern "C" {
// For configuration scripts to detect libais and the version numbers.
int LibAisVersionMajor();
int LibAisVersionMinor();
}

namespace libais {

// Returns vector of the text between the delimiters.  Uses an empty string
// for empty fields.  Empty string returns a vector of length 1 containing an
// empty string.
// Empty delim_str is not allowed.
vector<string> Split(const string &str, const string &delim_str);

// Returns the text in the nth field starting with the first field being 0.
// Empty delim_str is not allowed.
string GetNthField(const string &str, const size_t n, const string &delim_str);

// Returns the number of pad bits in an AIS AIVDM NMEA string.
// Returns -1 if there is an error.
int GetPad(const string &nmea_str);

// Returns the armored payload of an AIS AIVDM NMEA string.
// Returns an empty string if there was an error.
string GetBody(const string &nmea_str);


// Note: Needs to be kept in sync with AIS_STATUS_STRINGS list in ais.cpp.
enum AIS_STATUS {
  AIS_UNINITIALIZED,  // Message is not yet completely decoded.
  AIS_OK,
  AIS_ERR_BAD_BIT_COUNT,
  AIS_ERR_BAD_NMEA_CHR,
  AIS_ERR_BAD_PTR,
  AIS_ERR_UNKNOWN_MSG_TYPE,
  AIS_ERR_MSG_NOT_IMPLEMENTED,
  AIS_ERR_MSG_SUB_NOT_IMPLEMENTED,
  AIS_ERR_EXPECTED_STRING,
  AIS_ERR_BAD_MSG_CONTENT,
  AIS_ERR_MSG_TOO_LONG,
  AIS_ERR_BAD_SUB_MSG,
  AIS_ERR_BAD_SUB_SUB_MSG,
  AIS_STATUS_NUM_CODES
};

extern const char *const AIS_STATUS_STRINGS[AIS_STATUS_NUM_CODES];

// Designated Area Codes (DAC) / Maritime Identification Digits define
// which country controls a subset of the submessage spaces within
// AIS "binary" messages 6, 8, 25, and 26.  See:
//   http://www.itu.int/online/mms/glad/cga_mids.sh?lng=E
//   http://en.wikipedia.org/w/index.php?title=Maritime_identification_digits
// River Information System (RIS):
//   ECE-TRANS-SC3-2006-10r-RIS.pdf

enum Dac {
  AIS_DAC_0_TEST = 0,
  AIS_DAC_1_INTERNATIONAL,
  AIS_DAC_200_RIS = 200,
  AIS_DAC_201_ALBANIA = 201,
  AIS_DAC_202_ANDORRA = 202,
  AIS_DAC_203_AUSTRIA = 203,
  AIS_DAC_204_AZORES = 204,
  AIS_DAC_205_BELGIUM = 205,
  AIS_DAC_206_BELARUS = 206,
  AIS_DAC_207_BULGARIA = 207,
  AIS_DAC_208_VATICAN_CITY_STATE = 208,
  AIS_DAC_209_CYPRUS = 209,
  AIS_DAC_210_CYPRUS = 210,
  AIS_DAC_211_GERMANY = 211,
  AIS_DAC_212_CYPRUS = 212,
  AIS_DAC_213_GEORGIA = 213,
  AIS_DAC_214_MOLDOVA = 214,
  AIS_DAC_215_MALTA = 215,
  AIS_DAC_216_ARMENIA = 216,
  AIS_DAC_218_GERMANY = 218,
  AIS_DAC_219_DENMARK = 219,
  AIS_DAC_220_DENMARK = 220,
  AIS_DAC_224_SPAIN = 224,
  AIS_DAC_225_SPAIN = 225,
  AIS_DAC_226_FRANCE = 226,
  AIS_DAC_227_FRANCE = 227,
  AIS_DAC_228_FRANCE = 228,
  AIS_DAC_229_MALTA = 229,
  AIS_DAC_230_FINLAND = 230,
  AIS_DAC_231_FAROE_ISLANDS = 231,
  AIS_DAC_232_UNITED_KINGDOM_OF_GREAT_BRITAIN_AND_NORTHERN_IRELAND = 232,
  AIS_DAC_233_UNITED_KINGDOM_OF_GREAT_BRITAIN_AND_NORTHERN_IRELAND = 233,
  AIS_DAC_234_UNITED_KINGDOM_OF_GREAT_BRITAIN_AND_NORTHERN_IRELAND = 234,
  AIS_DAC_235_UNITED_KINGDOM_OF_GREAT_BRITAIN_AND_NORTHERN_IRELAND = 235,
  AIS_DAC_236_GIBRALTAR = 236,
  AIS_DAC_237_GREECE = 237,
  AIS_DAC_238_CROATIA = 238,
  AIS_DAC_239_GREECE = 239,
  AIS_DAC_240_GREECE = 240,
  AIS_DAC_241_GREECE = 241,
  AIS_DAC_242_MOROCCO = 242,
  AIS_DAC_243_HUNGARY = 243,
  AIS_DAC_244_NETHERLANDS = 244,
  AIS_DAC_245_NETHERLANDS = 245,
  AIS_DAC_246_NETHERLANDS = 246,
  AIS_DAC_247_ITALY = 247,
  AIS_DAC_248_MALTA = 248,
  AIS_DAC_249_MALTA = 249,
  AIS_DAC_250_IRELAND = 250,
  AIS_DAC_251_ICELAND = 251,
  AIS_DAC_252_LIECHTENSTEIN = 252,
  AIS_DAC_253_LUXEMBOURG = 253,
  AIS_DAC_254_MONACO = 254,
  AIS_DAC_255_MADEIRA = 255,
  AIS_DAC_256_MALTA = 256,
  AIS_DAC_257_NORWAY = 257,
  AIS_DAC_258_NORWAY = 258,
  AIS_DAC_259_NORWAY = 259,
  AIS_DAC_261_POLAND = 261,
  AIS_DAC_262_MONTENEGRO = 262,
  AIS_DAC_263_PORTUGAL = 263,
  AIS_DAC_264_ROMANIA = 264,
  AIS_DAC_265_SWEDEN = 265,
  AIS_DAC_266_SWEDEN = 266,
  AIS_DAC_267_SLOVAKIA = 267,
  AIS_DAC_268_SAN_MARINO = 268,
  AIS_DAC_269_SWITZERLAND = 269,
  AIS_DAC_270_CZECH_REPUBLIC = 270,
  AIS_DAC_271_TURKEY = 271,
  AIS_DAC_272_UKRAINE = 272,
  AIS_DAC_273_RUSSIAN_FEDERATION = 273,
  AIS_DAC_274_MACEDONIA = 274,
  AIS_DAC_275_LATVIA = 275,
  AIS_DAC_276_ESTONIA = 276,
  AIS_DAC_277_LITHUANIA = 277,
  AIS_DAC_278_SLOVENIA = 278,
  AIS_DAC_279_SERBIA = 279,
  AIS_DAC_301_ANGUILLA = 301,
  AIS_DAC_303_ALASKA = 303,
  AIS_DAC_304_ANTIGUA_AND_BARBUDA = 304,
  AIS_DAC_305_ANTIGUA_AND_BARBUDA = 305,
  AIS_DAC_306_NETHERLANDS_ANTILLES = 306,
  AIS_DAC_307_ARUBA = 307,
  AIS_DAC_308_BAHAMAS = 308,
  AIS_DAC_309_BAHAMAS = 309,
  AIS_DAC_310_BERMUDA = 310,
  AIS_DAC_311_BAHAMAS = 311,
  AIS_DAC_312_BELIZE = 312,
  AIS_DAC_314_BARBADOS = 314,
  AIS_DAC_316_CANADA = 316,
  AIS_DAC_319_CAYMAN_ISLANDS = 319,
  AIS_DAC_321_COSTA_RICA = 321,
  AIS_DAC_323_CUBA = 323,
  AIS_DAC_325_DOMINICA = 325,
  AIS_DAC_327_DOMINICAN_REPUBLIC = 327,
  AIS_DAC_329_GUADELOUPE = 329,
  AIS_DAC_330_GRENADA = 330,
  AIS_DAC_331_GREENLAND = 331,
  AIS_DAC_332_GUATEMALA = 332,
  AIS_DAC_334_HONDURAS = 334,
  AIS_DAC_336_HAITI = 336,
  AIS_DAC_338_UNITED_STATES_OF_AMERICA = 338,
  AIS_DAC_339_JAMAICA = 339,
  AIS_DAC_341_SAINT_KITTS_AND_NEVIS = 341,
  AIS_DAC_343_SAINT_LUCIA = 343,
  AIS_DAC_345_MEXICO = 345,
  AIS_DAC_347_MARTINIQUE = 347,
  AIS_DAC_348_MONTSERRAT = 348,
  AIS_DAC_350_NICARAGUA = 350,
  AIS_DAC_351_PANAMA = 351,
  AIS_DAC_352_PANAMA = 352,
  AIS_DAC_353_PANAMA = 353,
  AIS_DAC_354_PANAMA = 354,
  AIS_DAC_355_PANAMA = 355,
  AIS_DAC_356_PANAMA = 356,
  AIS_DAC_357_PANAMA = 357,
  AIS_DAC_358_PUERTO_RICO = 358,
  AIS_DAC_359_EL_SALVADOR = 359,
  AIS_DAC_361_SAINT_PIERRE_AND_MIQUELON = 361,
  AIS_DAC_362_TRINIDAD_AND_TOBAGO = 362,
  AIS_DAC_364_TURKS_AND_CAICOS_ISLANDS = 364,
  AIS_DAC_366_UNITED_STATES_OF_AMERICA = 366,
  AIS_DAC_367_UNITED_STATES_OF_AMERICA = 367,
  AIS_DAC_368_UNITED_STATES_OF_AMERICA = 368,
  AIS_DAC_369_UNITED_STATES_OF_AMERICA = 369,
  AIS_DAC_370_PANAMA = 370,
  AIS_DAC_371_PANAMA = 371,
  AIS_DAC_372_PANAMA = 372,
  AIS_DAC_373_PANAMA = 373,
  AIS_DAC_375_SAINT_VINCENT_AND_THE_GRENADINES = 375,
  AIS_DAC_376_SAINT_VINCENT_AND_THE_GRENADINES = 376,
  AIS_DAC_377_SAINT_VINCENT_AND_THE_GRENADINES = 377,
  AIS_DAC_378_BRITISH_VIRGIN_ISLANDS = 378,
  AIS_DAC_379_UNITED_STATES_VIRGIN_ISLANDS = 379,
  AIS_DAC_401_AFGHANISTAN = 401,
  AIS_DAC_403_SAUDI_ARABIA = 403,
  AIS_DAC_405_BANGLADESH = 405,
  AIS_DAC_408_BAHRAIN = 408,
  AIS_DAC_410_BHUTAN = 410,
  AIS_DAC_412_CHINA = 412,
  AIS_DAC_413_CHINA = 413,
  AIS_DAC_414_CHINA = 414,
  AIS_DAC_416_TAIWAN = 416,
  AIS_DAC_417_SRI_LANKA = 417,
  AIS_DAC_419_INDIA = 419,
  AIS_DAC_422_IRAN = 422,
  AIS_DAC_423_AZERBAIJAN = 423,
  AIS_DAC_425_IRAQ = 425,
  AIS_DAC_428_ISRAEL = 428,
  AIS_DAC_431_JAPAN = 431,
  AIS_DAC_432_JAPAN = 432,
  AIS_DAC_434_TURKMENISTAN = 434,
  AIS_DAC_436_KAZAKHSTAN = 436,
  AIS_DAC_437_UZBEKISTAN = 437,
  AIS_DAC_438_JORDAN = 438,
  AIS_DAC_440_KOREA = 440,
  AIS_DAC_441_KOREA = 441,
  AIS_DAC_443_STATE_OF_PALESTINE = 443,
  AIS_DAC_445_DEMOCRATIC_PEOPLES_REPUBLIC_OF_KOREA = 445,
  AIS_DAC_447_KUWAIT = 447,
  AIS_DAC_450_LEBANON = 450,
  AIS_DAC_451_KYRGYZ_REPUBLIC = 451,
  AIS_DAC_453_MACAO = 453,
  AIS_DAC_455_MALDIVES = 455,
  AIS_DAC_457_MONGOLIA = 457,
  AIS_DAC_459_NEPAL = 459,
  AIS_DAC_461_OMAN = 461,
  AIS_DAC_463_PAKISTAN = 463,
  AIS_DAC_466_QATAR = 466,
  AIS_DAC_468_SYRIAN_ARAB_REPUBLIC = 468,
  AIS_DAC_470_UNITED_ARAB_EMIRATES = 470,
  AIS_DAC_472_TAJIKISTAN = 472,
  AIS_DAC_473_YEMEN = 473,
  AIS_DAC_475_YEMEN = 475,
  AIS_DAC_477_HONG_KONG = 477,
  AIS_DAC_478_BOSNIA_AND_HERZEGOVINA = 478,
  AIS_DAC_501_ADELIE_LAND = 501,
  AIS_DAC_503_AUSTRALIA = 503,
  AIS_DAC_506_MYANMAR = 506,
  AIS_DAC_508_BRUNEI_DARUSSALAM = 508,
  AIS_DAC_510_MICRONESIA = 510,
  AIS_DAC_511_PALAU = 511,
  AIS_DAC_512_NEW_ZEALAND = 512,
  AIS_DAC_514_CAMBODIA = 514,
  AIS_DAC_515_CAMBODIA = 515,
  AIS_DAC_516_CHRISTMAS_ISLAND = 516,
  AIS_DAC_518_COOK_ISLANDS = 518,
  AIS_DAC_520_FIJI = 520,
  AIS_DAC_523_COCOS = 523,
  AIS_DAC_525_INDONESIA = 525,
  AIS_DAC_529_KIRIBATI = 529,
  AIS_DAC_531_LAO_PEOPLES_DEMOCRATIC_REPUBLIC = 531,
  AIS_DAC_533_MALAYSIA = 533,
  AIS_DAC_536_NORTHERN_MARIANA_ISLANDS = 536,
  AIS_DAC_538_MARSHALL_ISLANDS = 538,
  AIS_DAC_540_NEW_CALEDONIA = 540,
  AIS_DAC_542_NIUE = 542,
  AIS_DAC_544_NAURU = 544,
  AIS_DAC_546_FRENCH_POLYNESIA = 546,
  AIS_DAC_548_PHILIPPINES = 548,
  AIS_DAC_553_PAPUA_NEW_GUINEA = 553,
  AIS_DAC_555_PITCAIRN_ISLAND = 555,
  AIS_DAC_557_SOLOMON_ISLANDS = 557,
  AIS_DAC_559_AMERICAN_SAMOA = 559,
  AIS_DAC_561_SAMOA = 561,
  AIS_DAC_563_SINGAPORE = 563,
  AIS_DAC_564_SINGAPORE = 564,
  AIS_DAC_565_SINGAPORE = 565,
  AIS_DAC_566_SINGAPORE = 566,
  AIS_DAC_567_THAILAND = 567,
  AIS_DAC_570_TONGA = 570,
  AIS_DAC_572_TUVALU = 572,
  AIS_DAC_574_VIETNAM = 574,
  AIS_DAC_576_VANUATU = 576,
  AIS_DAC_577_VANUATU = 577,
  AIS_DAC_578_WALLIS_AND_FUTUNA_ISLANDS = 578,
  AIS_DAC_601_SOUTH_AFRICA = 601,
  AIS_DAC_603_ANGOLA = 603,
  AIS_DAC_605_ALGERIA = 605,
  AIS_DAC_607_SAINT_PAUL_AND_AMSTERDAM_ISLANDS = 607,
  AIS_DAC_608_ASCENSION_ISLAND = 608,
  AIS_DAC_609_BURUNDI = 609,
  AIS_DAC_610_BENIN = 610,
  AIS_DAC_611_BOTSWANA = 611,
  AIS_DAC_612_CENTRAL_AFRICAN_REPUBLIC = 612,
  AIS_DAC_613_CAMEROON = 613,
  AIS_DAC_615_CONGO = 615,
  AIS_DAC_616_COMOROS = 616,
  AIS_DAC_617_CABO_VERDE = 617,
  AIS_DAC_618_CROZET_ARCHIPELAGO = 618,
  AIS_DAC_619_COTE_DIVOIRE = 619,
  AIS_DAC_620_COMOROS = 620,
  AIS_DAC_621_DJIBOUTI = 621,
  AIS_DAC_622_EGYPT = 622,
  AIS_DAC_624_ETHIOPIA = 624,
  AIS_DAC_625_ERITREA = 625,
  AIS_DAC_626_GABONESE_REPUBLIC = 626,
  AIS_DAC_627_GHANA = 627,
  AIS_DAC_629_GAMBIA = 629,
  AIS_DAC_630_GUINEABISSAU = 630,
  AIS_DAC_631_EQUATORIAL_GUINEA = 631,
  AIS_DAC_632_GUINEA = 632,
  AIS_DAC_633_BURKINA_FASO = 633,
  AIS_DAC_634_KENYA = 634,
  AIS_DAC_635_KERGUELEN_ISLANDS = 635,
  AIS_DAC_636_LIBERIA = 636,
  AIS_DAC_637_LIBERIA = 637,
  AIS_DAC_638_SOUTH_SUDAN = 638,
  AIS_DAC_642_LIBYA = 642,
  AIS_DAC_644_LESOTHO = 644,
  AIS_DAC_645_MAURITIUS = 645,
  AIS_DAC_647_MADAGASCAR = 647,
  AIS_DAC_649_MALI = 649,
  AIS_DAC_650_MOZAMBIQUE = 650,
  AIS_DAC_654_MAURITANIA = 654,
  AIS_DAC_655_MALAWI = 655,
  AIS_DAC_656_NIGER = 656,
  AIS_DAC_657_NIGERIA = 657,
  AIS_DAC_659_NAMIBIA = 659,
  AIS_DAC_660_REUNION = 660,
  AIS_DAC_661_RWANDA = 661,
  AIS_DAC_662_SUDAN = 662,
  AIS_DAC_663_SENEGAL = 663,
  AIS_DAC_664_SEYCHELLES = 664,
  AIS_DAC_665_SAINT_HELENA = 665,
  AIS_DAC_666_SOMALIA = 666,
  AIS_DAC_667_SIERRA_LEONE = 667,
  AIS_DAC_668_SAO_TOME_AND_PRINCIPE = 668,
  AIS_DAC_669_SWAZILAND = 669,
  AIS_DAC_670_CHAD = 670,
  AIS_DAC_671_TOGOLESE_REPUBLIC = 671,
  AIS_DAC_672_TUNISIA = 672,
  AIS_DAC_674_TANZANIA = 674,
  AIS_DAC_675_UGANDA = 675,
  AIS_DAC_676_DEMOCRATIC_REPUBLIC_OF_THE_CONGO = 676,
  AIS_DAC_677_TANZANIA = 677,
  AIS_DAC_678_ZAMBIA = 678,
  AIS_DAC_679_ZIMBABWE = 679,
  AIS_DAC_701_ARGENTINE_REPUBLIC = 701,
  AIS_DAC_710_BRAZIL = 710,
  AIS_DAC_720_BOLIVIA = 720,
  AIS_DAC_725_CHILE = 725,
  AIS_DAC_730_COLOMBIA = 730,
  AIS_DAC_735_ECUADOR = 735,
  AIS_DAC_740_FALKLAND_ISLANDS = 740,
  AIS_DAC_745_GUIANA = 745,
  AIS_DAC_750_GUYANA = 750,
  AIS_DAC_755_PARAGUAY = 755,
  AIS_DAC_760_PERU = 760,
  AIS_DAC_765_SURINAME = 765,
  AIS_DAC_770_URUGUAY = 770,
  AIS_DAC_775_VENEZUELA = 775
};

class AisBitset;

class AisPoint {
 public:
  double lng_deg;
  double lat_deg;

  AisPoint();
  AisPoint(double lng_deg_, double lat_deg_);
};
ostream& operator<< (ostream &o, const AisPoint &position);

class AisMsg {
 public:
  int message_id;
  int repeat_indicator;
  int mmsi;

  // TODO(schwehr): make status private and have accessors.
  bool had_error() const {  return status != AIS_OK;  }
  AIS_STATUS get_error() const { return status; }

  virtual ~AisMsg() {}

 protected:
  AIS_STATUS status;  // AIS_OK or error code
  int num_chars;  // Number of characters in the nmea_payload.
  size_t num_bits;  // Number of bits in the nmea_payload.

  AisMsg() : status(AIS_UNINITIALIZED), num_chars(0), num_bits(0) {}
  AisMsg(const char *nmea_payload, const size_t pad);
};

// TODO(schwehr): factor out commstate from all messages?
class Ais1_2_3 : public AisMsg {
 public:
  int nav_status;
  bool rot_over_range;
  int rot_raw;
  float rot;
  float sog;  // Knots.
  int position_accuracy;
  AisPoint position;
  float cog;  // Degrees.
  int true_heading;
  // TODO(schwehr): What about a leap second when timestamp 60 may be valid?
  int timestamp;
  int special_manoeuvre;
  int spare;
  bool raim;

  // COMM state SOTDMA msgs 1 and 2
  int sync_state;  // SOTDMA and ITDMA
  bool slot_timeout_valid;
  int slot_timeout;

  // Based on slot_timeout which ones are valid
  bool received_stations_valid;
  int received_stations;

  bool slot_number_valid;
  int slot_number;

  bool utc_valid;
  int utc_hour;
  int utc_min;
  int utc_spare;

  bool slot_offset_valid;
  int slot_offset;

  // ITDMA - msg type 3
  bool slot_increment_valid;
  int slot_increment;

  bool slots_to_allocate_valid;
  int slots_to_allocate;

  bool keep_flag_valid;
  bool keep_flag;  // 3.3.7.3.2 Annex 2 ITDMA.  Table 20

  Ais1_2_3(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais1_2_3 &msg);

// 4 bsreport and 11 utc date response
class Ais4_11 : public AisMsg {
 public:
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  int position_accuracy;
  AisPoint position;
  int fix_type;
  int transmission_ctl;
  int spare;
  bool raim;

  // COMM state SOTDMA msgs 1 and 2
  int sync_state;
  int slot_timeout;

  // Based on slot_timeout which ones are valid
  bool received_stations_valid;
  int received_stations;

  bool slot_number_valid;
  int slot_number;

  bool utc_valid;
  int utc_hour;
  int utc_min;
  int utc_spare;

  bool slot_offset_valid;
  int slot_offset;

  // **NO** ITDMA
  Ais4_11(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais4_11 &msg);

class Ais5 : public AisMsg {
 public:
  int ais_version;
  int imo_num;
  string callsign;
  string name;
  int type_and_cargo;
  int dim_a;
  int dim_b;
  int dim_c;
  int dim_d;
  int fix_type;
  int eta_month;
  int eta_day;
  int eta_hour;
  int eta_minute;
  float draught;  // present static draft. m
  string destination;
  int dte;
  int spare;

  Ais5(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais5 &msg);

// AIS Binary Broadcast message ... parent to many
class Ais6 : public AisMsg {
 public:
  int seq;  // sequence number
  int mmsi_dest;
  bool retransmit;
  int spare;
  int dac;  // dac+fi = app id
  int fi;

  // TODO(schwehr): how to make Ais6 protected?
  Ais6(const char *nmea_payload, const size_t pad);

 protected:
  Ais6() {}
};
ostream& operator<< (ostream &o, const Ais6 &msg);

// http://www.e-navigation.nl/content/monitoring-aids-navigation
// Zeni Lite Buoy Co., Ltd buoy status.
class Ais6_0_0 : public Ais6 {
 public:
  int sub_id;
  float voltage;
  float current;
  bool dc_power_supply;  // False is AC.
  bool light_on;
  bool battery_low;
  bool off_position;
  int spare2;

  Ais6_0_0(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_0_0 &msg);

// Text message.  ITU 1371-1
class Ais6_1_0 : public Ais6 {
 public:
  bool ack_required;
  int msg_seq;
  string text;
  int spare2;

  Ais6_1_0(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_0 &msg);

// Application ack.  ITU 1371-1
class Ais6_1_1 : public Ais6 {
 public:
  int ack_dac;
  int msg_seq;
  int spare2;

  Ais6_1_1(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_1 &msg);

// Interrogation for a DAC/FI.  ITU 1371-1
class Ais6_1_2 : public Ais6 {
 public:
  int req_dac;
  int req_fi;
  // TODO(schwehr): spare2?

  Ais6_1_2(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_2 &msg);

// Capability interogation.  ITU 1371-1
class Ais6_1_3 : public Ais6 {
 public:
  int req_dac;
  int spare2;

  Ais6_1_3(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_3 &msg);

// Capability interogation reply.  ITU 1371-1
class Ais6_1_4 : public Ais6 {
 public:
  int ack_dac;
  std::array<int, 64> capabilities;
  std::array<int, 64> cap_reserved;
  int spare2;

  Ais6_1_4(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_4 &msg);

// International function message 5: Application ack to addr binary message.
class Ais6_1_5 : public Ais6 {
 public:
  // TODO(schwehr): How to handle the sequence number and retransmit flag?
  int ack_dac;
  int ack_fi;
  int seq_num;
  bool ai_available;  // TODO(schwehr): AI?  Is this the dac/fi being acked?
  int ai_response;
  int spare;

  Ais6_1_5(const char *nmea_payload, const size_t pad);
};

// IMO Circ 236 Dangerous cargo indication
// Not to be transmitted after 2012-Jan-01
class Ais6_1_12 : public Ais6 {
 public:
  string last_port;
  int utc_month_dep;  // actual time of departure
  int utc_day_dep;
  int utc_hour_dep;
  int utc_min_dep;
  string next_port;
  int utc_month_next;  // estimated arrival
  int utc_day_next;
  int utc_hour_next;
  int utc_min_next;
  string main_danger;
  string imo_cat;
  int un;
  int value;  // TODO(schwehr): units?
  int value_unit;
  int spare2;

  Ais6_1_12(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_12 &msg);


class Ais6_1_14_Window {
 public:
  AisPoint position;
  int utc_hour_from;
  int utc_min_from;
  int utc_hour_to;
  int utc_min_to;
  int cur_dir;
  float cur_speed;
};


// IMO Circ 236 Tidal window
// Not to be transmitted after 2012-Jan-01
class Ais6_1_14 : public Ais6 {
 public:
  int utc_month;
  int utc_day;
  vector<Ais6_1_14_Window> windows;

  Ais6_1_14(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_14 &msg);


// IMO Circ 289 Clearance time to enter port
class Ais6_1_18 : public Ais6 {
 public:
  int link_id;
  int utc_month;
  int utc_day;
  int utc_hour;
  int utc_min;
  string port_berth;
  string dest;
  AisPoint position;
  std::array<int, 2> spare2;  // 32 bits per spare

  Ais6_1_18(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_18 &msg);


// IMO Circ 289 Berthing data
class Ais6_1_20 : public Ais6 {
 public:
  int link_id;
  int length;
  float depth;
  int mooring_position;
  int utc_month;
  int utc_day;
  int utc_hour;
  int utc_min;
  bool services_known;
  // TODO(schwehr): enum of service types
  std::array<int, 26> services;
  string name;
  AisPoint position;

  Ais6_1_20(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_20 &msg);

class Ais6_1_25_Cargo {
 public:
  int code_type;
  bool imdg_valid;  // also set with BC
  int imdg;
  bool spare_valid;
  int spare;  // imdg or dc or marpols
  bool un_valid;
  int un;
  bool bc_valid;
  int bc;
  bool marpol_oil_valid;
  int marpol_oil;
  bool marpol_cat_valid;
  int marpol_cat;

  Ais6_1_25_Cargo();
  // TODO(schwehr): Add a constructor from an AisBitset.
};

// IMO Circ 289 Dangerous cargo indication 2
// Replaces 8_1_12?
class Ais6_1_25 : public Ais6 {
 public:
  int amount_unit;
  int amount;

  vector<Ais6_1_25_Cargo> cargos;  // 0 to 17 cargo entries

  Ais6_1_25(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_25 &msg);

// TODO(schwehr): Addressed sensor report 6_1_26.
// TODO(schwehr): IMO Circ 289 Route information 6_1_28.
// TODO(schwehr): IMO Circ 289 Text description 6_1_30.

// IMO Circ 289
// Warning: The bit encoding for 6_1_14_Window and 6_1_32 on
//   the wire has x and y in a different order.
// TODO(schwehr): Reuse Ais6_1_14_Window
class Ais6_1_32_Window {
 public:
  AisPoint position;
  int from_utc_hour;
  int from_utc_min;
  int to_utc_hour;
  int to_utc_min;
  int cur_dir;
  float cur_speed;  // knots

  Ais6_1_32_Window();
};


// IMO Circ 289 Tidal window
class Ais6_1_32 : public Ais6 {
 public:
  int utc_month;
  int utc_day;
  vector<Ais6_1_32_Window> windows;

  Ais6_1_32(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_32 &msg);

// Number of persons on board.  ITU 1371-1
class Ais6_1_40 : public Ais6 {
 public:
  int persons;
  int spare2;

  Ais6_1_40(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais6_1_40 &msg);

//////////////////////////////////////////////////////////////////////

// 7 and 13 are ACKs for msg 6 and 12
class Ais7_13 : public AisMsg {
 public:
  int spare;

  vector<int> dest_mmsi;
  vector<int> seq_num;

  Ais7_13(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais7_13 &msg);

// AIS Binary Broadcast message ... parent to many
class Ais8 : public AisMsg {
 public:
  int spare;
  // TODO(schwehr): seq? // ITU M.R. 1371-3 Anex 2 5.3.1
  int dac;  // dac+fi = app id
  int fi;

  // TODO(schwehr): make Ais8 protected
  Ais8(const char *nmea_payload, const size_t pad);

 protected:
  Ais8() {}
};
ostream& operator<< (ostream &o, const Ais8 &msg);

// Text telegram ITU 1371-1
class Ais8_1_0 : public Ais8 {
 public:
  bool ack_required;
  int msg_seq;
  string text;
  int spare2;

  Ais8_1_0(const char *nmea_payload, size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_0 &msg);

// 8_1_1 No message
// 8_1_2 No message
// 8_1_3 No message
// 8_1_4 No message

// IMO Circ 289 met hydro - Not to be transmitted after 2013-Jan-01
// See also IMO Circ 236
class Ais8_1_11 : public Ais8 {
 public:
  AisPoint position;
  int day;
  int hour;
  int minute;
  int wind_ave;  // kts
  int wind_gust;  // kts
  int wind_dir;
  int wind_gust_dir;
  float air_temp;  // C
  int rel_humid;
  float dew_point;
  float air_pres;
  int air_pres_trend;
  float horz_vis;  // NM
  float water_level;  // m
  int water_level_trend;
  float surf_cur_speed;
  int surf_cur_dir;
  float cur_speed_2;  // kts
  int cur_dir_2;
  int cur_depth_2;  // m
  float cur_speed_3;  // kts
  int cur_dir_3;
  int cur_depth_3;  // m
  float wave_height;  // m
  int wave_period;
  int wave_dir;
  float swell_height;
  int swell_period;
  int swell_dir;
  int sea_state;  // beaufort scale
  float water_temp;
  int precip_type;
  float salinity;  // Part per mil (1/1000).
  int ice;  // yes/no/undef/unknown
  int spare2;
  int extended_water_level;  // OHMEX uses this for extra water level precision

  Ais8_1_11(const char *nmea_payload, size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_11 &msg);


// IMO Circ 236 Fairway closed - Not to be transmitted after 2012-Jan-01
class Ais8_1_13 : public Ais8 {
 public:
  string reason;
  string location_from;
  string location_to;
  int radius;
  int units;
  // TODO(schwehr): utc?  warning: day/month out of order
  int day_from;
  int month_from;
  int hour_from;
  int minute_from;

  int day_to;
  int month_to;
  int hour_to;
  int minute_to;
  int spare2;

  Ais8_1_13(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_13 &msg);

// IMO Circ 236 Extended ship static and voyage data
// Not to be transmitted after 2012-Jan-01
class Ais8_1_15 : public Ais8 {
 public:
  float air_draught;
  int spare2;

  Ais8_1_15(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_15 &msg);

// IMO Circ 236 Number of persons on board
class Ais8_1_16 : public Ais8 {
 public:
  int persons;
  int spare2;

  Ais8_1_16(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_16 &msg);

class Ais8_1_17_Target {
 public:
  int type;
  string id;
  int spare;
  AisPoint position;
  int cog;
  int timestamp;
  int sog;

  Ais8_1_17_Target();
};

// IMO Circ 236 VTS Generated/synthetic targets
class Ais8_1_17 : public Ais8 {
 public:
  vector<Ais8_1_17_Target> targets;

  Ais8_1_17(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_17 &msg);

// No 8_1_18

// IMO Circ 289 Marine traffic signal
class Ais8_1_19 : public Ais8 {
 public:
  int link_id;
  string name;
  AisPoint position;  // funny bit count
  int status;
  int signal;
  int utc_hour_next;
  int utc_min_next;
  int next_signal;
  std::array<int, 4> spare2;

  Ais8_1_19(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_19 &msg);

// No message 8_1_20

// IMO Circ 289 Weather observation report from ship
class Ais8_1_21 : public Ais8 {
 public:
  int type_wx_report;

  // TYPE 0
  string location;
  AisPoint position;  // 25, 24 bits
  int utc_day;
  int utc_hour;
  int utc_min;
  // wx - use wx[0]
  float horz_viz;  // nautical miles
  int humidity;  // %
  int wind_speed;  // ave knots
  int wind_dir;
  float pressure;  // hPa - float needed for type 1
  int pressure_tendency;
  float air_temp;  // C
  float water_temp;  // C
  int wave_period;  // s
  float wave_height;
  int wave_dir;
  float swell_height;  // m
  int swell_dir;
  int swell_period;  // s
  int spare2;

  // TYPE 1 - !@#$!!!!!
  // x, y
  int utc_month;
  // utc_day, hour, min
  int cog;
  float sog;
  int heading;  // Assume this is true degrees????
  // pressure defined in type 0
  float rel_pressure;  // 3 hour hPa
  // pressure_tendenc defined in type 0
  // wind_dir defined in type 0
  float wind_speed_ms;  // m/s
  int wind_dir_rel;
  float wind_speed_rel;  // m/s
  float wind_gust_speed;  // m/s
  int wind_gust_dir;
  int air_temp_raw;  // TODO(schwehr): Convert this to C.  Kelvin makes no sense
  // humidity defined in type 0
  // sea_temp_k
  int water_temp_raw;  // TODO(schwehr): fix this
  // hor_viz
  int wx[3];  // current, past 1, past 2
  int cloud_total;
  int cloud_low;
  int cloud_low_type;
  int cloud_middle_type;
  int cloud_high_type;
  float alt_lowest_cloud_base;
  // wave_period
  // wave_height
  // swell_dir
  // swell_period
  // swell_height
  int swell_dir_2;
  int swell_period_2;
  float swell_height_2;
  float ice_thickness;  // Network is cm, storing m.
  int ice_accretion;
  int ice_accretion_cause;
  int sea_ice_concentration;
  int amt_type_ice;
  int ice_situation;
  int ice_devel;
  int bearing_ice_edge;

  Ais8_1_21(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_21 &msg);


const size_t AIS8_1_22_NUM_NAMES = 128;
const size_t AIS8_1_22_SUBAREA_SIZE = 87;
extern const char *ais8_1_22_notice_names[AIS8_1_22_NUM_NAMES];

enum Ais8_1_22_AreaShapeEnum {
  AIS8_1_22_SHAPE_ERROR = -1,
  AIS8_1_22_SHAPE_CIRCLE = 0,  // OR Point.
  AIS8_1_22_SHAPE_RECT = 1,
  AIS8_1_22_SHAPE_SECTOR = 2,
  AIS8_1_22_SHAPE_POLYLINE = 3,
  AIS8_1_22_SHAPE_POLYGON = 4,
  AIS8_1_22_SHAPE_TEXT = 5,
  AIS8_1_22_SHAPE_RESERVED_6 = 6,
  AIS8_1_22_SHAPE_RESERVED_7 = 7
};

extern const char *ais8_1_22_shape_names[8];

// Sub-Areas for the Area Notice class

class Ais8_1_22_SubArea {
 public:
  virtual Ais8_1_22_AreaShapeEnum getType() const = 0;
  virtual ~Ais8_1_22_SubArea() {}
};

Ais8_1_22_SubArea*
ais8_1_22_subarea_factory(const AisBitset &bs,
                          const size_t offset);

// or Point if radius is 0
class Ais8_1_22_Circle : public Ais8_1_22_SubArea {
 public:
  AisPoint position;  // Longitude and latitude.
  // Going to assume that the precision is not useful.
  int precision;  // How many decimal places for x and y.
  int radius_m;
  unsigned int spare;  // 18 bits.

  Ais8_1_22_Circle(const AisBitset &bs, const size_t offset);
  ~Ais8_1_22_Circle() {}
  Ais8_1_22_AreaShapeEnum getType() const {return AIS8_1_22_SHAPE_CIRCLE;}
};

class Ais8_1_22_Rect : public Ais8_1_22_SubArea {
 public:
  AisPoint position;  // Longitude and latitude.
  int precision;  // How many decimal places for x and y.  Useless.
  int e_dim_m;  // East dimension in meters.
  int n_dim_m;
  int orient_deg;  // Orientation in degrees from true north.
  unsigned int spare;  // 5 bits.

  Ais8_1_22_Rect(const AisBitset &bs, const size_t offset);
  ~Ais8_1_22_Rect() {}
  Ais8_1_22_AreaShapeEnum getType() const {return AIS8_1_22_SHAPE_RECT;}
};

class Ais8_1_22_Sector : public Ais8_1_22_SubArea {
 public:
  AisPoint position;  // Longitude and latitude.
  // TODO(schwehr): precision in IMO, but not RTCM.  Double check.
  int precision;  // How many decimal places for x and y?
  int radius_m;
  int left_bound_deg;
  int right_bound_deg;

  Ais8_1_22_Sector(const AisBitset &bs, const size_t offset);
  ~Ais8_1_22_Sector() {}
  Ais8_1_22_AreaShapeEnum getType() const {return AIS8_1_22_SHAPE_SECTOR;}
};

// Or Waypoint
// Must have a point before on the VDL
// TODO(schwehr): do I bring in the prior point x, y, precision?
class Ais8_1_22_Polyline : public Ais8_1_22_SubArea {
 public:
  // TODO(schwehr): int precision; // How many decimal places for x and y.

  // Up to 4 points
  vector<float> angles;
  vector<float> dists_m;
  unsigned int spare;  // 2 bit.

  Ais8_1_22_Polyline(const AisBitset &bs, const size_t offset);
  ~Ais8_1_22_Polyline() {}
  Ais8_1_22_AreaShapeEnum getType() const {return AIS8_1_22_SHAPE_POLYLINE;}
};

// TODO(schwehr): Bring in the prior point?  And do we fold the sub area data
// into one polygon if there are more than one?
class Ais8_1_22_Polygon : public Ais8_1_22_SubArea {
 public:
  // TODO(schwehr): int precision; // How many decimal places for x and y.

  // Up to 4 points in a first message, but aggregated if multiple sub areas
  vector<float> angles;
  vector<float> dists_m;
  unsigned int spare;  // 2 bit

  Ais8_1_22_Polygon(const AisBitset &bs, const size_t offset);
  ~Ais8_1_22_Polygon() {}
  Ais8_1_22_AreaShapeEnum getType() const {return AIS8_1_22_SHAPE_POLYGON;}
};


class Ais8_1_22_Text : public Ais8_1_22_SubArea {
 public:
  string text;
  // TODO(schwehr): spare?

  Ais8_1_22_Text(const AisBitset &bs, const size_t offset);
  ~Ais8_1_22_Text() {}
  Ais8_1_22_AreaShapeEnum getType() const {return AIS8_1_22_SHAPE_TEXT;}
};

// Area Notice class

class Ais8_1_22 : public Ais8 {
 public:
  int link_id;  // 10 bit id to match up text blocks.
  int notice_type;  // Area type / Notice Description.
  int month;  // These are in UTC.
  int day;
  int hour;
  int minute;
  int duration_minutes;  // Time from the start until the notice expires.

  // 1 or more sub messages
  vector<Ais8_1_22_SubArea *> sub_areas;

  Ais8_1_22(const char *nmea_payload, const size_t pad);
  ~Ais8_1_22();
};
ostream& operator<< (ostream& o, Ais8_1_22 const& msg);

// No message 8_1_23

// IMO Circ 289 Extended ship static and voyage-related
class Ais8_1_24 : public Ais8 {
 public:
  int link_id;
  float air_draught;  // m
  string last_port;
  std::array<std::string, 2> next_ports;

  // TODO(schwehr): enum list of param types
  std::array<int, 26> solas_status;  // 0 NA, 1 operational, 2 SNAFU, 3 no data
  int ice_class;
  int shaft_power;  // horses
  int vhf;
  string lloyds_ship_type;
  int gross_tonnage;
  int laden_ballast;
  int heavy_oil;
  int light_oil;
  int diesel;
  int bunker_oil;  // tonnes
  int persons;
  int spare2;

  Ais8_1_24(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_24 &msg);

// No message 8_1_25

const size_t AIS8_1_26_REPORT_SIZE = 112;

enum Ais8_1_26_SensorEnum {
  AIS8_1_26_SENSOR_ERROR = -1,
  AIS8_1_26_SENSOR_LOCATION = 0,
  AIS8_1_26_SENSOR_STATION = 1,
  AIS8_1_26_SENSOR_WIND = 2,
  AIS8_1_26_SENSOR_WATER_LEVEL = 3,
  AIS8_1_26_SENSOR_CURR_2D = 4,
  AIS8_1_26_SENSOR_CURR_3D = 5,
  AIS8_1_26_SENSOR_HORZ_FLOW = 6,
  AIS8_1_26_SENSOR_SEA_STATE = 7,
  AIS8_1_26_SENSOR_SALINITY = 8,
  AIS8_1_26_SENSOR_WX = 9,
  AIS8_1_26_SENSOR_AIR_DRAUGHT = 10,
  AIS8_1_26_SENSOR_RESERVED_11 = 11,
  AIS8_1_26_SENSOR_RESERVED_12 = 12,
  AIS8_1_26_SENSOR_RESERVED_13 = 13,
  AIS8_1_26_SENSOR_RESERVED_14 = 14,
  AIS8_1_26_SENSOR_RESERVED_15 = 15,
};

class Ais8_1_26_SensorReport {
 public:
  int report_type;
  int utc_day;
  int utc_hr;
  int utc_min;
  int site_id;  // aka link_id

  virtual Ais8_1_26_SensorEnum getType() const = 0;
  virtual ~Ais8_1_26_SensorReport() {}
};

Ais8_1_26_SensorReport*
ais8_1_26_sensor_report_factory(const AisBitset &bs,
                                const size_t offset);

class Ais8_1_26_Location : public Ais8_1_26_SensorReport {
 public:
  AisPoint position;
  float z;  // alt in m from MSL
  int owner;
  int timeout;
  int spare;

  Ais8_1_26_Location(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Location() {}
  Ais8_1_26_SensorEnum getType() const { return AIS8_1_26_SENSOR_LOCATION; }
};

class Ais8_1_26_Station : public Ais8_1_26_SensorReport {
 public:
  string name;
  int spare;

  Ais8_1_26_Station(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Station() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_STATION;}
};

class Ais8_1_26_Wind : public Ais8_1_26_SensorReport {
 public:
  int wind_speed;  // knots
  int wind_gust;  // knots
  int wind_dir;
  int wind_gust_dir;
  int sensor_type;
  int wind_forecast;  // knots
  int wind_gust_forecast;  // knots
  int wind_dir_forecast;
  int utc_day_forecast;
  int utc_hour_forecast;
  int utc_min_forecast;
  int duration;
  int spare;

  Ais8_1_26_Wind(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Wind() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_WIND;}
};

class Ais8_1_26_WaterLevel : public Ais8_1_26_SensorReport {
 public:
  int type;
  float level;  // m.  assuming it is being stored at 0.01 m inc.
  int trend;
  int vdatum;
  int sensor_type;
  int forecast_type;
  float level_forecast;
  int utc_day_forecast;
  int utc_hour_forecast;
  int utc_min_forecast;
  int duration;  // minutes
  int spare;

  Ais8_1_26_WaterLevel(const AisBitset &bs, const size_t offset);
  Ais8_1_26_WaterLevel() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_WATER_LEVEL;}
};

class Ais8_1_26_Curr2D_Current {
 public:
  float speed;  // knots
  int dir;
  int depth;  // m
};

class Ais8_1_26_Curr2D : public Ais8_1_26_SensorReport {
 public:
  Ais8_1_26_Curr2D_Current currents[3];
  int type;
  int spare;

  Ais8_1_26_Curr2D(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Curr2D() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_CURR_2D;}
};

class Ais8_1_26_Curr3D_Current {
 public:
  float north;
  float east;
  float up;
  int depth;  // m
};

class Ais8_1_26_Curr3D : public Ais8_1_26_SensorReport {
 public:
  Ais8_1_26_Curr3D_Current currents[2];
  int type;
  int spare;

  Ais8_1_26_Curr3D(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Curr3D() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_CURR_3D;}
};

class Ais8_1_26_HorzFlow_Current {
 public:
  int bearing;  // deg
  int dist;  // m
  float speed;  // knots
  int dir;  // deg
  int level;  // m
};

class Ais8_1_26_HorzFlow : public Ais8_1_26_SensorReport {
 public:
  Ais8_1_26_HorzFlow_Current currents[2];
  int spare;

  Ais8_1_26_HorzFlow(const AisBitset &bs, const size_t offset);
  Ais8_1_26_HorzFlow() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_HORZ_FLOW;}
};

class Ais8_1_26_SeaState : public Ais8_1_26_SensorReport {
 public:
  float swell_height;
  int swell_period;  // seconds
  int swell_dir;  // deg
  int sea_state;
  int swell_sensor_type;
  float water_temp;  // C
  float water_temp_depth;  // m
  int water_sensor_type;
  float wave_height;
  int wave_period;  // seconds
  int wave_dir;  // deg
  int wave_sensor_type;
  float salinity;

  Ais8_1_26_SeaState(const AisBitset &bs, const size_t offset);
  Ais8_1_26_SeaState() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_SEA_STATE;}
};

class Ais8_1_26_Salinity : public Ais8_1_26_SensorReport {
 public:
  float water_temp;  // C
  float conductivity;  // siemens/m
  float pressure;  // decibars
  float salinity;  // 0/00 ppt
  int salinity_type;
  int sensor_type;
  int spare[2];

  Ais8_1_26_Salinity(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Salinity() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_SALINITY;}
};

class Ais8_1_26_Wx : public Ais8_1_26_SensorReport {
 public:
  float air_temp;  // C
  int air_temp_sensor_type;
  int precip;
  float horz_vis;  // nm
  float dew_point;  // C
  int dew_point_type;
  float air_pressure;  // Pascals (Pa).
  int air_pressure_trend;
  int air_pressor_type;
  float salinity;  // 0/00 ppt
  int spare;

  Ais8_1_26_Wx(const AisBitset &bs, const size_t offset);
  Ais8_1_26_Wx() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_WX;}
};

class Ais8_1_26_AirDraught : public Ais8_1_26_SensorReport {
 public:
  float draught;
  float gap;
  float forecast_gap;
  int trend;
  int utc_day_forecast;
  int utc_hour_forecast;
  int utc_min_forecast;
  int spare;

  Ais8_1_26_AirDraught(const AisBitset &bs, const size_t offset);
  Ais8_1_26_AirDraught() {}
  Ais8_1_26_SensorEnum getType() const {return AIS8_1_26_SENSOR_AIR_DRAUGHT;}
};

// IMO Circ 289 Environmental
class Ais8_1_26 : public Ais8 {
 public:
  vector<Ais8_1_26_SensorReport *> reports;

  Ais8_1_26(const char *nmea_payload, const size_t pad);
  ~Ais8_1_26();
};
ostream& operator<< (ostream &o, const Ais8_1_26 &msg);


// IMO Circ 289 Route information
class Ais8_1_27 : public Ais8 {
 public:
  int link_id;
  int sender_type;
  int route_type;
  int utc_month;
  int utc_day;
  int utc_hour;
  int utc_min;
  int duration;
  vector<AisPoint> waypoints;

  Ais8_1_27(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_27 &msg);


//  No message 8_1_28


// IMO Circ 289 Text description
class Ais8_1_29 : public Ais8 {
 public:
  int link_id;
  string text;
  int spare2;

  Ais8_1_29(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_29 &msg);


// No message 8_1_30

// IMO Circ 289 Meteorological and Hydrographic data
// Section 1, Table 1.1
// TODO(schwehr): is this exactly the same as 8_1_11 or has anything changed?
//       x,y swapped.
class Ais8_1_31 : public Ais8 {
 public:
  AisPoint position;  // Opposite the bit order of 8_1_11
  int position_accuracy;  // New field
  int utc_day;
  int utc_hour;
  int utc_min;
  int wind_ave;  // kts
  int wind_gust;  // kts
  int wind_dir;
  int wind_gust_dir;
  float air_temp;  // C
  int rel_humid;
  float dew_point;
  float air_pres;  // Pascals (Pa).
  int air_pres_trend;
  float horz_vis;  // NM
  float water_level;  // m
  int water_level_trend;

  float surf_cur_speed;
  int surf_cur_dir;
  float cur_speed_2;  // kts
  int cur_dir_2;
  int cur_depth_2;  // m
  float cur_speed_3;  // kts
  int cur_dir_3;
  int cur_depth_3;  // m

  float wave_height;  // m
  int wave_period;
  int wave_dir;
  float swell_height;  // m
  int swell_period;
  int swell_dir;
  int sea_state;  // beaufort scale - Table 1.2
  float water_temp;
  int precip_type;
  float salinity;
  int ice;  // yes/no/undef/unknown
  int spare2;

  Ais8_1_31(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais8_1_31 &msg);

// TODO(schwehr): Persons on board ITU 1371-1 8_1_40.

// ECE-TRANS-SC3-2006-10e-RIS.pdf - River Information System
// Inland ship static and voyage related data
class Ais8_200_10 : public Ais8 {
 public:
  string eu_id;  // European Vessel ID - 8 characters
  float length;  // m
  float beam;  // m
  int ship_type;
  int haz_cargo;
  float draught;
  int loaded;
  // Sensor quality.
  int speed_qual;
  int course_qual;
  int heading_qual;
  int spare2;

  Ais8_200_10(const char *nmea_payload, const size_t pad);
};

// 21 and 22 do not exist

// ECE-TRANS-SC3-2006-10e-RIS.pdf - River Information System
class Ais8_200_23 : public Ais8 {
 public:
  int utc_year_start;
  int utc_month_start;
  int utc_day_start;
  int utc_year_end;
  int utc_month_end;
  int utc_day_end;
  int utc_hour_start;
  int utc_min_start;
  int utc_hour_end;
  int utc_min_end;
  AisPoint position1;
  AisPoint position2;
  int type;
  int min;
  int max;
  int classification;
  int wind_dir;  // EMMA CODE
  int spare2;

  Ais8_200_23(const char *nmea_payload, const size_t pad);
};

// ECE-TRANS-SC3-2006-10e-RIS.pdf - River Information System
// Water Level
class Ais8_200_24 : public Ais8 {
 public:
  string country;
  std::array<int, 4> gauge_ids;
  std::array<float, 4> levels;  // m

  Ais8_200_24(const char *nmea_payload, const size_t pad);
};

// ECE-TRANS-SC3-2006-10e-RIS.pdf - River Information System
class Ais8_200_40 : public Ais8 {
 public:
  AisPoint position;
  int form;
  int dir;  // degrees
  int stream_dir;
  int status_raw;
  // TODO(schwehr): int status[9];  // WTF is the encoding for this?
  int spare2;

  Ais8_200_40(const char *nmea_payload, const size_t pad);
};

// ECE-TRANS-SC3-2006-10e-RIS.pdf - River Information System
class Ais8_200_55 : public Ais8 {
 public:
  int crew;
  int passengers;
  int yet_more_personnel;
  std::array<int, 3> spare2;  // 51 spare bits.

  Ais8_200_55(const char *nmea_payload, const size_t pad);
};

enum Ais8_366_22_AreaShapeEnum {
  AIS8_366_22_SHAPE_ERROR = -1,
  AIS8_366_22_SHAPE_CIRCLE = 0,
  AIS8_366_22_SHAPE_RECT = 1,
  AIS8_366_22_SHAPE_SECTOR = 2,
  AIS8_366_22_SHAPE_POLYLINE = 3,
  AIS8_366_22_SHAPE_POLYGON = 4,
  AIS8_366_22_SHAPE_TEXT = 5,
  AIS8_366_22_SHAPE_RESERVED_6 = 6,
  AIS8_366_22_SHAPE_RESERVED_7 = 7
};

extern const char *shape_names[8];

class Ais8_366_22_SubArea {
 public:
    virtual Ais8_366_22_AreaShapeEnum getType() = 0;
    virtual ~Ais8_366_22_SubArea() { }
};

Ais8_366_22_SubArea*
ais8_366_22_subarea_factory(const AisBitset &bs,
                            const size_t offset);

// or Point if radius is 0
class Ais8_366_22_Circle : public Ais8_366_22_SubArea {
 public:
  AisPoint position;
  // TODO(schwehr): int precision
  int radius_m;
  unsigned int spare;

  Ais8_366_22_Circle(const AisBitset &bs, const size_t offset);
  ~Ais8_366_22_Circle() {}
  Ais8_366_22_AreaShapeEnum getType() {return AIS8_366_22_SHAPE_CIRCLE;}
};

class Ais8_366_22_Rect : public Ais8_366_22_SubArea {
 public:
  AisPoint position;  // longitude and latitude
  // TODO(schwehr): int precision
  int e_dim_m;  // East dimension in meters
  int n_dim_m;
  int orient_deg;  // Orientation in degrees from true north
  unsigned int spare;  // 5 bits

  Ais8_366_22_Rect(const AisBitset &bs, const size_t offset);
  ~Ais8_366_22_Rect() {}
  Ais8_366_22_AreaShapeEnum getType() {return AIS8_366_22_SHAPE_RECT;}
};

class Ais8_366_22_Sector : public Ais8_366_22_SubArea {
 public:
  AisPoint position;
  // TODO(schwehr): int precision
  int radius_m;
  int left_bound_deg;
  int right_bound_deg;
  // TODO(schwehr): spare?

  Ais8_366_22_Sector(const AisBitset &bs, const size_t offset);
  ~Ais8_366_22_Sector() {}
  Ais8_366_22_AreaShapeEnum getType() {return AIS8_366_22_SHAPE_SECTOR;}
};

// Or Waypoint
// Must have a point before on the VDL, but pulled together here.
class Ais8_366_22_Polyline : public Ais8_366_22_SubArea {
 public:
  AisPoint position;  // longitude and latitude
  // TODO(schwehr): precision

  // Up to 4 points
  vector<float> angles;
  vector<float> dists_m;
  unsigned int spare;

  Ais8_366_22_Polyline(const AisBitset &bs, const size_t offset);
  ~Ais8_366_22_Polyline() {}
  Ais8_366_22_AreaShapeEnum getType() {return AIS8_366_22_SHAPE_POLYLINE;}
};

class Ais8_366_22_Polygon : public Ais8_366_22_SubArea {
 public:
  AisPoint position;  // longitude and latitude
  // TODO(schwehr): precision?

  // Up to 4 points in a first message, but aggregated if multiple sub areas
  vector<float> angles;
  vector<float> dists_m;
  unsigned int spare;

  Ais8_366_22_Polygon(const AisBitset &bs, const size_t offset);
  ~Ais8_366_22_Polygon() {}
  Ais8_366_22_AreaShapeEnum getType() {return AIS8_366_22_SHAPE_POLYGON;}
};

class Ais8_366_22_Text : public Ais8_366_22_SubArea {
 public:
  string text;
  unsigned int spare;  // 3 bits

  Ais8_366_22_Text(const AisBitset &bs, const size_t offset);
  ~Ais8_366_22_Text() {}
  Ais8_366_22_AreaShapeEnum getType() {return AIS8_366_22_SHAPE_TEXT;}
};

class Ais8_366_22 : public Ais8 {
 public:
  // Common block at the front
  int link_id;  // 10 bit id to match up text blocks
  int notice_type;  // area_type / Notice Description
  int month;  // These really are in utc
  int day;
  int utc_hour;
  int utc_minute;
  int duration_minutes;  // Time from the start until the notice expires
  // 1 or more sub messages

  vector<Ais8_366_22_SubArea *> sub_areas;

  Ais8_366_22(const char *nmea_payload, const size_t pad);
  ~Ais8_366_22();
};
ostream& operator<< (ostream& o, Ais8_366_22 const& msg);

const size_t AIS8_366_22_NUM_NAMES = 128;
extern const char *ais8_366_22_notice_names[AIS8_366_22_NUM_NAMES];

// 366 34 - Kurt older whale message 2008-2010
// TODO(schwehr): Ais8_366_34

// USCG Blue Force encrypted message.
class Ais8_366_56 : public Ais8 {
 public:
  vector<unsigned char> encrypted;

  Ais8_366_56(const char *nmea_payload, const size_t pad);
};

class Ais8_367_22_SubArea {
 public:
  virtual Ais8_366_22_AreaShapeEnum getType() const = 0;
  virtual ~Ais8_367_22_SubArea() { }
};

Ais8_367_22_SubArea*
ais8_367_22_subarea_factory(const AisBitset &bs,
                            const size_t offset);

class Ais8_367_22_Circle : public Ais8_367_22_SubArea {
 public:
  AisPoint position;
  int precision;
  int radius_m;
  unsigned int spare;

  Ais8_367_22_Circle(const AisBitset &bs, const size_t offset);
  ~Ais8_367_22_Circle() {}
  Ais8_366_22_AreaShapeEnum getType() const {return AIS8_366_22_SHAPE_CIRCLE;}
};

class Ais8_367_22_Rect : public Ais8_367_22_SubArea {
 public:
  AisPoint position;
  int precision;
  int e_dim_m;
  int n_dim_m;
  int orient_deg;
  unsigned int spare;

  Ais8_367_22_Rect(const AisBitset &bs, const size_t offset);
  ~Ais8_367_22_Rect() {}
  Ais8_366_22_AreaShapeEnum getType() const {return AIS8_366_22_SHAPE_RECT;}
};

class Ais8_367_22_Sector : public Ais8_367_22_SubArea {
 public:
  AisPoint position;
  int precision;
  int radius_m;
  int left_bound_deg;
  int right_bound_deg;
  int spare;

  Ais8_367_22_Sector(const AisBitset &bs, const size_t offset);
  ~Ais8_367_22_Sector() {}
  Ais8_366_22_AreaShapeEnum getType() const {return AIS8_366_22_SHAPE_SECTOR;}
};

// Polyline or Polygon
class Ais8_367_22_Poly : public Ais8_367_22_SubArea {
 public:
  Ais8_366_22_AreaShapeEnum shape;
  AisPoint position;
  int precision;

  // Up to 4 points
  vector<float> angles;
  vector<float> dists_m;
  unsigned int spare;

  Ais8_367_22_Poly(const AisBitset &bs, const size_t offset,
                   Ais8_366_22_AreaShapeEnum area_shape);
  ~Ais8_367_22_Poly() {}
  Ais8_366_22_AreaShapeEnum getType() const {return shape;}
};

class Ais8_367_22_Text : public Ais8_367_22_SubArea {
 public:
  string text;
  unsigned int spare;  // 3 bits

  Ais8_367_22_Text(const AisBitset &bs, const size_t offset);
  ~Ais8_367_22_Text() {}
  Ais8_366_22_AreaShapeEnum getType() const {return AIS8_366_22_SHAPE_TEXT;}
};

class Ais8_367_22 : public Ais8 {
 public:
  int version;
  int link_id;
  int notice_type;
  int month;
  int day;
  int hour;
  int minute;
  int duration_minutes;
  int spare2;

  vector<Ais8_367_22_SubArea *> sub_areas;

  Ais8_367_22(const char *nmea_payload, const size_t pad);
  ~Ais8_367_22();
};
ostream& operator<< (ostream& o, Ais8_367_22 const& msg);

class Ais9 : public AisMsg {
 public:
  int alt;  // m above sea level
  int sog;
  int position_accuracy;
  AisPoint position;
  float cog;
  int timestamp;
  int alt_sensor;
  int spare;
  int dte;
  int spare2;
  int assigned_mode;
  bool raim;
  int commstate_flag;

  int sync_state;  // In both SOTDMA and ITDMA

  // SOTDMA
  bool slot_timeout_valid;
  int slot_timeout;

  // Based on slot_timeout which ones are valid
  bool received_stations_valid;
  int received_stations;

  bool slot_number_valid;
  int slot_number;

  bool utc_valid;
  int utc_hour;
  int utc_min;
  int utc_spare;

  bool slot_offset_valid;
  int slot_offset;

  // ITDMA
  bool slot_increment_valid;
  int slot_increment;

  bool slots_to_allocate_valid;
  int slots_to_allocate;

  bool keep_flag_valid;
  bool keep_flag;

  Ais9(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais9 &msg);

// 10 ":" - UTC and date inquiry
class Ais10 : public AisMsg {
 public:
  int spare;
  int dest_mmsi;
  int spare2;

  Ais10(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais10 &msg);

// 11 ';' - See 4_11

// 12 - '<' - Addressed safety related text.
class Ais12 : public AisMsg {
 public:
  int seq_num;
  int dest_mmsi;
  bool retransmitted;
  int spare;
  string text;
  int spare2;

  Ais12(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais12 &msg);

// 13 '=' - See 7

// 14 - '>' - Safety broadcast text.
class Ais14 : public AisMsg {
 public:
  int spare;
  string text;
  int spare2;

  Ais14(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais14 &msg);

// 15 - '?' - Interrogation
class Ais15 : public AisMsg {
 public:
  int spare;
  int mmsi_1;
  int msg_1_1;
  int slot_offset_1_1;

  int spare2;
  int dest_msg_1_2;
  int slot_offset_1_2;

  int spare3;
  int mmsi_2;
  int msg_2;
  int slot_offset_2;
  int spare4;

  Ais15(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais15 &msg);

// 16 - '@' - Assigned mode command
class Ais16 : public AisMsg {
 public:
  int spare;
  int dest_mmsi_a;
  int offset_a;
  int inc_a;
  int dest_mmsi_b;
  int offset_b;
  int inc_b;
  int spare2;

  Ais16(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais16 &msg);

// ITU-R M.823  http://www.itu.int/rec/R-REC-M.823/en
// 17 - 'A' - GNSS broacast - TODO(schwehr): only partially coded
class Ais17 : public AisMsg {
 public:
  int spare;
  AisPoint position;
  int spare2;
  int gnss_type;
  int z_cnt;
  int station;
  int seq;
  // N - do not need to store this
  int health;
  // TODO(schwehr): Handle payload

  Ais17(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais17 &msg);


// 18 - 'B' - Class B position report
class Ais18 : public AisMsg {
 public:
  int spare;
  float sog;
  int position_accuracy;
  AisPoint position;  // Long and lat
  float cog;
  int true_heading;
  int timestamp;
  int spare2;
  int unit_flag;  // 0 is SOTDMA, 1 is Carrier Sense (CS).
  int display_flag;
  int dsc_flag;
  int band_flag;
  int m22_flag;
  int mode_flag;
  bool raim;
  int commstate_flag;

  // SOTDMA
  int sync_state;
  bool slot_timeout_valid;
  int slot_timeout;

  // Based on slot_timeout which ones are valid
  bool received_stations_valid;
  int received_stations;

  bool slot_number_valid;
  int slot_number;

  bool utc_valid;
  int utc_hour;
  int utc_min;
  int utc_spare;

  bool slot_offset_valid;
  int slot_offset;

  // ITDMA
  bool slot_increment_valid;
  int slot_increment;

  bool slots_to_allocate_valid;
  int slots_to_allocate;

  bool keep_flag_valid;
  bool keep_flag;

  // If commstate set to 1 for Carrier Sense (CS) devices, there is no
  // state and the commstate region is supposed to be filled with this value:
  //   1100000000000000110
  bool commstate_cs_fill_valid;
  int commstate_cs_fill;

  Ais18(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais18 &msg);

// 19 - 'C' - Class B extended ship and position
class Ais19 : public AisMsg {
 public:
  int spare;
  float sog;
  int position_accuracy;
  AisPoint position;  // Long and lat
  float cog;
  int true_heading;
  int timestamp;
  int spare2;
  string name;
  int type_and_cargo;
  int dim_a;
  int dim_b;
  int dim_c;
  int dim_d;
  int fix_type;
  bool raim;
  int dte;
  int assigned_mode;
  int spare3;

  Ais19(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais19 &msg);

// 20 - 'D' - Data link management
// TODO(schwehr): consider a vector
class Ais20 : public AisMsg {
 public:
  int spare;
  int offset_1;
  int num_slots_1;
  int timeout_1;
  int incr_1;

  bool group_valid_2;
  int offset_2;
  int num_slots_2;
  int timeout_2;
  int incr_2;

  bool group_valid_3;
  int offset_3;
  int num_slots_3;
  int timeout_3;
  int incr_3;

  bool group_valid_4;
  int offset_4;
  int num_slots_4;
  int timeout_4;
  int incr_4;

  int spare2;

  Ais20(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais20 &msg);

// 21 - 'E' - Aids to navigation report
class Ais21 : public AisMsg {
 public:
  int aton_type;
  string name;
  int position_accuracy;
  AisPoint position;
  int dim_a;
  int dim_b;
  int dim_c;
  int dim_d;
  int fix_type;
  int timestamp;
  bool off_pos;
  int aton_status;
  bool raim;
  bool virtual_aton;
  bool assigned_mode;
  int spare;
  // Extended name goes on the end of name
  int spare2;

  Ais21(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais21 &msg);

// 22 - 'F' - Channel Management
class Ais22 : public AisMsg {
 public:
  int spare;
  int chan_a;
  int chan_b;
  int txrx_mode;
  bool power_low;

  // if addressed false, then geographic position
  bool pos_valid;
  AisPoint position1;
  AisPoint position2;

  // if addressed is true
  bool dest_valid;
  int dest_mmsi_1;
  int dest_mmsi_2;

  int chan_a_bandwidth;
  int chan_b_bandwidth;
  int zone_size;

  int spare2;  // Lame that they make a huge spare here.  Bad bad bad

  Ais22(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais22 &msg);

// 23 - 'G' - Group Assignment Command
class Ais23 : public AisMsg {
 public:
  int spare;
  AisPoint position1;
  AisPoint position2;
  int station_type;
  int type_and_cargo;

  int spare2;  // 22 bits of spare here?  what were people thinking?

  int txrx_mode;
  int interval_raw;  // raw value, not sec
  int quiet;
  int spare3;

  Ais23(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais23 &msg);

// 24 - 'H' - Class B Static Data report
class Ais24 : public AisMsg {
 public:
  int part_num;

  // Part A
  string name;

  // Part B
  int type_and_cargo;
  string vendor_id;
  string callsign;
  int dim_a;
  int dim_b;
  int dim_c;
  int dim_d;
  int spare;

  // Part C - Not defined by ITU 1371-5
  // Part D - Not defined by ITU 1371-5

  Ais24(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais24 &msg);

// 25 - 'I' - Single slot binary message - addressed or broadcast
// TODO(schwehr): handle payload
class Ais25 : public AisMsg {
 public:
  bool use_app_id;  // if false, payload is unstructured binary.

  bool dest_mmsi_valid;
  int dest_mmsi;  // only valid if addressed
  // If unstructured:
  // TODO(schwehr): vector<unsigned char> payload;

  int dac;  // valid if use_app_id is true
  int fi;

  Ais25(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais25 &msg);

// 26 - 'J' - Multi slot binary message with comm state
// TODO(schwehr): handle payload
class Ais26 : public AisMsg {
 public:
  bool use_app_id;  // if false, payload is unstructured binary

  bool dest_mmsi_valid;
  int dest_mmsi;  // only valid if addressed

  int dac;  // valid it use_app_id
  int fi;

  // TODO(schwehr): vector<unsigned char> payload;  // If unstructured.  Yuck.

  int commstate_flag;  // 0 - SOTDMA, 1 - TDMA

  // SOTDMA
  int sync_state;
  bool slot_timeout_valid;
  int slot_timeout;

  // Based on slot_timeout which ones are valid
  bool received_stations_valid;
  int received_stations;

  bool slot_number_valid;
  int slot_number;

  bool utc_valid;
  int utc_hour;
  int utc_min;
  int utc_spare;

  bool slot_offset_valid;
  int slot_offset;

  // ITDMA
  bool slot_increment_valid;
  int slot_increment;

  bool slots_to_allocate_valid;
  int slots_to_allocate;

  bool keep_flag_valid;
  bool keep_flag;

  Ais26(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais26 &msg);

// 27 - 'K' - Long-range position report - e.g. for satellite receivers
class Ais27 : public AisMsg {
 public:
  int position_accuracy;
  bool raim;
  int nav_status;
  AisPoint position;
  int sog;  // Knots.
  int cog;  // Degrees.
  bool gnss;  // warning: bits in AIS are flipped sense
  int spare;

  Ais27(const char *nmea_payload, const size_t pad);
};
ostream& operator<< (ostream &o, const Ais27 &msg);

//////////////////////////////////////////////////////////////////////
// Support class for decoding
//////////////////////////////////////////////////////////////////////
static const int MAX_BITS = 1192;

class AisBitset : protected bitset<MAX_BITS> {
 public:
  AisBitset();

  AIS_STATUS ParseNmeaPayload(const char *nmea_payload, int pad);

  int GetNumBits() const { return num_bits; }
  int GetNumChars() const { return num_chars; }
  int GetPosition() const { return current_position; }
  int GetRemaining() const { return num_bits - current_position; }

  const AisBitset& SeekRelative(int d) const;
  const AisBitset& SeekTo(size_t pos) const;

  bool operator[](size_t pos) const;

  unsigned int ToUnsignedInt(const size_t start, const size_t len) const;
  int ToInt(const size_t start, const size_t len) const;
  string ToString(const size_t start, const size_t len) const;

  const AisPoint ToAisPoint(const size_t start, const size_t point_size) const;

  // Visible for testing.
  static bitset<6> Reverse(const bitset<6> &bits);

 protected:
  // TODO(schwehr): do not use long
  typedef union {
    long long_val;  // NOLINT
    unsigned long ulong_val;  // NOLINT
  } long_union;

  int num_bits;
  int num_chars;

  static bool nmea_ord_initialized_;
  static bitset<6> nmea_ord_[128];
  static const char bits_to_char_tbl_[];

  static void InitNmeaOrd();

 private:
  // This will help uncover dicontinuities when querying sequential bits, i.e.
  // when we query a bit sequence that is not in direct succession of the
  // previous one. In the future, we may use this to automatically determine
  // the next read location.
  // This field is also used to determine the number of remaining bits after the
  // last read position.
  mutable int current_position;
};


}  // namespace libais

#endif  // LIBAIS_AIS_H_
