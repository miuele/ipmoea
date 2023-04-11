#ifndef IPMOEA_COMPOSITOR_H_INCLUDED
#define IPMOEA_COMPOSITOR_H_INCLUDED

#include <utility>

namespace ipmoea {

namespace utils {

template <class T, class F, class ...Funcs>
struct compositor_impl : compositor_impl<T, Funcs...> {
	constexpr compositor_impl(F &&f, Funcs &&...funcs)
		: compositor_impl<T, Funcs...>(std::forward<Funcs>(funcs)...), f_(std::forward<F>(f))
	{
	}

	constexpr T operator()(T value) {
		return f_(compositor_impl<T, Funcs...>::operator()(value));
	}

	F f_;
};

template <class T, class F>
struct compositor_impl<T, F> {
	constexpr compositor_impl(F &&f)
		: f_(std::forward<F>(f))
	{}

	constexpr T operator()(T value) {
		return f_(value);
	}

	F f_;
};

template <class T, class ...Funcs>
struct compositor
	: private compositor_impl<T, Funcs...>
{
	constexpr compositor(Funcs &&...funcs)
		: compositor_impl<T, Funcs...>(std::forward<Funcs>(funcs)...)
	{
	}

	constexpr T operator()(T value) {
		return compositor_impl<T, Funcs...>::operator()(value);
	}
};

template <class T, class ...Funcs>
constexpr auto compose(Funcs &&...funcs) {
	return compositor<T, Funcs...>(std::forward<Funcs>(funcs)...);
}

}

}

#endif
