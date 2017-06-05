#include "altruct/io/iostream_overloads.h"

#include <cmath>
#include <stdint.h>

iostream_fraction_state_t iostream_fraction_state;
bool iostream_fraction_state_init =
    io_fraction_denominator::set_default() &&
    io_fraction_as_pair::set_default();

iostream_modulo_state_t iostream_modulo_state;
bool iostream_modulo_state_init =
    io_modulo_modulus::set_default() &&
    io_modulo_as_pair::set_default();

iostream_polynom_state_t iostream_polynom_state;
bool iostream_polynom_state_init =
    io_polynom_as_vector::set_default();
