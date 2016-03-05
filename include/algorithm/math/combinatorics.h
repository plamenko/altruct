#pragma once

#include <algorithm>
#include <utility>

namespace altruct {
namespace math {

/**
 * Rearranges the elements in the range [begin, mid) and [mid, end) into the next lexicographically greater combination.
 * Selected elements are in the [begin, mid), the rest are in the [mid, end).
 */
template<typename BidIt>
bool next_combination(BidIt begin, BidIt mid, BidIt end) {
	BidIt rd, wr;
	for (BidIt cur = mid; cur != begin;) {
		for (BidIt cmp = wr = --cur; wr != mid; cmp = wr, ++wr) {
			for (rd = wr; rd != end && !(*cmp < *rd); ++rd);
			if (rd != end) {
				std::swap(*wr, *rd);
				// bubble(wr, rd, end)
				for (BidIt next = rd; --next, next != wr && *rd < *next; rd = next) std::swap(*next, *rd);
				for (BidIt next = rd; ++next, next != end && *next < *rd; rd = next) std::swap(*next, *rd);
			} else {
				++cmp;
				// merge(cur, cmp, end)
				std::sort(cur, end);
				break;
			}
		}
		if (wr == mid) return true;
	}
	return false;
}

/**
 * Rearranges the elements in the range [begin, end) into the next lexicographically smaller partition.
 */
template<typename BidIt>
bool next_partition(BidIt begin, BidIt end) {
	typedef typename std::iterator_traits<BidIt>::value_type I;
	I l = 0, x = 0;
	for (BidIt it = end; it != begin;) {
		--it; l = l + 1;
 		x += *it;
		I xm = *it - 1;
		if (xm > 0 && (x - 1) / xm <= (l - 1)) {
			for (; x > xm; x -= xm) {
				*it++ = xm;
			}
			*it++ = x;
			return true;
		}
		*it = 0;
	}
	*begin = x;
	return false;
}

} // math
} // altruct
