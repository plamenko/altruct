HEADERS = subdir_glob([
    #('include', '**/*.h'),
    ('include', '**/base.h'),
    ('include', '**/modulo.h'),
  ],
  excludes = [
    'include/concurrency/concurrency.h',
  ],
  prefix='');

cxx_library(
  name = 'altruct_lib',
  srcs = glob([
    #'lib/**/*.cpp',
    'lib/**/base.cpp',
    'lib/**/modulo.cpp',
  ],
  excludes=[
  ]),
  headers = HEADERS,
)

cxx_test(
  name = 'altruct_tests',
  srcs = glob([
    #'test/**/*.cpp',
    'test/**/base_test.cpp',
    'test/**/modulo_test.cpp',
  ],
  excludes=[
    'test/concurrency/concurrency_test.cpp'
  ]),
  headers = HEADERS,
  deps = [
    ':altruct_lib',
  ],
)
