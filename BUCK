HEADERS = subdir_glob([
    ('include', '**/*.h'),
    ('experimental/include', '**/*.h'),
  ],
  excludes = [
    'include/concurrency/concurrency.h',
  ],
  prefix='');

TEST_UTIL_HEADERS = subdir_glob([
    ('test_util/include', '**/*.h'),
  ],
  prefix='');

def merge_dicts(list_of_dicts):
  r = {}
  for d in list_of_dicts:
    r.update(d)
  return r;

cxx_library(
  name = 'altruct_lib',
  srcs = glob([
    'lib/**/*.cpp',
  ],
  excludes=[
  ]),
  headers = HEADERS,
)

cxx_test(
  name = 'altruct_tests',
  srcs = glob([
    'test/**/*.cpp',
  ],
  excludes=[
    'test/concurrency/concurrency_test.cpp'
  ]),
  headers = merge_dicts([HEADERS, TEST_UTIL_HEADERS]),
  deps = [
    ':altruct_lib',
  ],
)

cxx_binary(
  name = 'altruct_samples',
  srcs = glob([
    'sample/**/*.cpp',
  ],
  excludes=[
  ]),
  headers = HEADERS,
  deps = [
    ':altruct_lib',
  ],
)
