#include "algorithm/random/xorshift.h"

#include "gtest/gtest.h"

using namespace std;
using namespace altruct::random;

uint64_t kTest1024_1_seed = 281878792946930618ULL;
uint64_t kTest1024_1_expected[] = {
	5680888935564682062ULL,
	12187295268616620767ULL,
	4505827850632960666ULL,
	5825356497907253918ULL,
	14260800475362949259ULL,
	12940355228812596651ULL,
	12228866189743246434ULL,
	14858741442510275433ULL,
	13766918008986895688ULL,
	11883179133175603530ULL,
	15893287646107542167ULL,
	15673802190103605763ULL,
	7141212171506088397ULL,
	72924962998261864ULL,
	10530881423460142722ULL,
	15639974282192402942ULL,
	9686174083177471174ULL,
	18090518010393227249ULL,
	7813995086651459506ULL
};

uint64_t kTest1024_2_seed[16] = {
	777894986665880512ULL,
	720101215243649686ULL,
	619717220485127537ULL,
	638523201128980145ULL,
	797654180308422169ULL,
	486366638561994011ULL,
	825831472224064496ULL,
	214730964401385301ULL,
	111622135404622164ULL,
	274965130298052407ULL,
	901667313599927378ULL,
	530022705104910973ULL,
	926477170849097080ULL,
	151756505878391145ULL,
	106798143325768204ULL,
	128775920927847452ULL
};
uint64_t kTest1024_2_expected[] = {
	2885546577937574243ULL,
	10495692899075972505ULL,
	17432112644418960426ULL,
	15166231084493567322ULL,
	15053372882255691589ULL,
	12533126124981720964ULL,
	3348845571680027832ULL,
	16824316703249453083ULL,
	12355388557269016206ULL,
	11128622746893912568ULL,
	2324202292288266896ULL,
	10507053014425423325ULL,
	6827273240345279207ULL,
	15525547251534011962ULL,
	13795729555081120273ULL,
	2846976465278434753ULL,
	4268598077660268637ULL,
	18004872390376912708ULL,
	7215581113654252882ULL
};

uint64_t kTest64_1_seed = 281878792946930618ULL;
uint64_t kTest64_1_expected[] = {
	7286020999113445271ULL,
	1679787891260431913ULL,
	17184147381894586086ULL,
	9215257752446926913ULL,
	5503049314823406013ULL,
	4485928818749946006ULL,
	5602597868101118508ULL,
	712687439541726861ULL,
	15148967293808197328ULL,
	4888420282158065575ULL,
	11027067857486436681ULL,
	11872143816678099104ULL,
	16048753931747004901ULL,
	1397531627844889849ULL,
	10770221634859117494ULL,
	11566780580534883098ULL,
	9680280807611523822ULL,
	14827137907540074284ULL,
	3986790697632744448ULL
};

uint64_t kTest64_2_seed = 901667313599927378ULL;
uint64_t kTest64_2_expected[] = {
	13801716981810879454ULL,
	2694637068706158203ULL,
	6666217419296894864ULL,
	11103292563332182035ULL,
	466352475997224035ULL,
	13574044771525795780ULL,
	13083901015177442753ULL,
	17399697188495373907ULL,
	13596084698471223298ULL,
	9726054091418748638ULL,
	10657446282345482264ULL,
	6346056096209814721ULL,
	5399111073589621412ULL,
	14503172178559587663ULL,
	11237479319194397685ULL,
	8936921881012674636ULL,
	1546044010778503679ULL,
	3883684463055549597ULL,
	16595650600837003373ULL
};

TEST(xorshift_test, default_constructor) {
	xorshift_1024star rng;
	for (int i = 0; i < 20; i++) {
		EXPECT_EQ(0, rng.next());
	}
}

TEST(xorshift_test, default_constructor_later_seed_short) {
	xorshift_1024star rng;
	rng.seed(kTest1024_1_seed);
	for (auto expected : kTest1024_1_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_test, default_constructor_later_seed_full) {
	xorshift_1024star rng;
	rng.seed(kTest1024_2_seed);
	for (auto expected : kTest1024_2_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_test, short_seed_constructor) {
	xorshift_1024star rng(kTest1024_1_seed);
	for (auto expected : kTest1024_1_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_test, full_seed_constructor) {
	xorshift_1024star rng(kTest1024_2_seed);
	for (auto expected : kTest1024_2_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_test, reseed) {
	xorshift_1024star rng(kTest1024_1_seed);
	rng.seed(kTest1024_2_seed);
	for (auto expected : kTest1024_2_expected) {
		EXPECT_EQ(expected, rng.next());
	}
	rng.seed(kTest1024_1_seed);
	for (auto expected : kTest1024_1_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_test, output_1000) {
	xorshift_1024star rng(kTest1024_1_seed);
	uint64_t v = 0;
	for (int i = 0; i < 1000; i++) {
		v = rng.next();
	}
	EXPECT_EQ(563674104727552105ULL, v);
}

TEST(xorshift_test, next_range) {
	xorshift_1024star rng(kTest1024_1_seed);
	EXPECT_EQ(162, rng.next(100, 1100 - 1));
	EXPECT_EQ(867, rng.next(100, 1100 - 1));
}

TEST(xorshift_test, next_uniform) {
	xorshift_1024star rng(kTest1024_1_seed);
	EXPECT_EQ(162, rng.next_uniform(100, 1100 - 1));
	EXPECT_EQ(867, rng.next_uniform(100, 1100 - 1));
}

TEST(xorshift_test, next_0_1) {
	xorshift_1024star rng(kTest1024_1_seed);
	EXPECT_NEAR(0.307961606279405, rng.next_0_1(), 1e-14);
}

TEST(xorshift_64star_test, default_constructor) {
	xorshift_64star rng;
	for (int i = 0; i < 20; i++) {
		EXPECT_EQ(0, rng.next());
	}
}

TEST(xorshift_64star_test, default_constructor_later_seed) {
	xorshift_64star rng;
	rng.seed(kTest64_1_seed);
	for (auto expected : kTest64_1_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_64star_test, short_seed_constructor) {
	xorshift_64star rng(kTest64_1_seed);
	for (auto expected : kTest64_1_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}

TEST(xorshift_64star_test, reseed) {
	xorshift_64star rng(kTest64_1_seed);
	rng.seed(kTest64_2_seed);
	for (auto expected : kTest64_2_expected) {
		EXPECT_EQ(expected, rng.next());
	}
	rng.seed(kTest64_1_seed);
	for (auto expected : kTest64_1_expected) {
		EXPECT_EQ(expected, rng.next());
	}
}
