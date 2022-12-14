/* https://learn.microsoft.com/en-us/cpp/error-messages/compiler-warnings/compiler-warnings-by-compiler-version?view=msvc-170
 */
#if defined(_MSC_VER)
#pragma warning(disable : 5240) // ignored attribute
#pragma warning(disable : 4649) // attributes are ignored in this context
//# if _MSC_VER > _MSC_SOME_VERSION
//#  pragma warning(disable: xxxx) // disable one more for special version
//# endif
#elif defined(__clang__)
#pragma clang diagnostic ignored "-Wignored-attributes"
//#  if __has_warning("-Wnew-special-warning")
//#   pragma clang diagnostic ignored "-Wnew-special-warning"
//#  endif
#elif defined(__GNUG__)
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
