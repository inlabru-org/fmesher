/* Idea from https://stackoverflow.com/questions/4193476/is-using-pragma-warning-push-pop-the-right-way-to-temporarily-alter-warning-lev */

#if defined(_MSC_VER)
# pragma warning(push)
#elif defined(__clang__)
# pragma clang diagnostic push
#elif defined(__GNUG__)
# pragma GCC diagnostic push
#endif
