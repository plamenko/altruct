HEADERS = subdir_glob([
    ('include', '**/*.h'),
  ],
  excludes = [
    'include/concurrency/concurrency.h',
  ],
  prefix='');

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
  headers = HEADERS,
  deps = [
    ':altruct_lib',
  ],
)
