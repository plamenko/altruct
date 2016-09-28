#pragma once

#define ALTRUCT_STRINGIFY(x) #x
#define ALTRUCT_TOSTRING(x) ALTRUCT_STRINGIFY(x)
#define ALTRUCT_AT __FILE__ ":" ALTRUCT_TOSTRING(__LINE__)
