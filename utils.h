#ifndef IPMOEA_UTILS_H_INCLUDED
#define IPMOEA_UTILS_H_INCLUDED

#include "zsequence.h"

namespace ipmoea {

namespace utils {

template <class T>
struct iir_result_set {
	T value;
};

template <class T>
struct iir_result_not_set {
	constexpr iir_result_set<T> operator=(T value) {
		return { value };
	}
};

template <class T, int N>
struct iir_input;

template <class T, int N>
struct iir_output;

namespace detail {

template <class T, int N, int I>
constexpr const T &iir_input_z(const iir_input<T, N> &o, std::integral_constant<int, I>) {
	return o.seq_.template top<(-I - 1)>();
}

template <class T, int N>
constexpr const T &iir_input_z(const iir_input<T, N> &o, std::integral_constant<int, 0>) {
	return o.top_;
}

template <class T, int N, int I>
constexpr const T &iir_output_z(const iir_output<T, N> &o, std::integral_constant<int, I>) {
	return o.seq_.template top<(-I - 1)>();
}

template <class T, int N>
constexpr iir_result_not_set<T> iir_output_z(const iir_output<T, N> &o, std::integral_constant<int, 0>) {
	return iir_result_not_set<T>{};
}

}

using ::zseq::zsequence;

template <class T, int N>
struct iir_input {

	constexpr explicit iir_input(T top, const zsequence<T, N> &past_input)
		: top_(top), seq_(past_input)
	{}

	template <int I>
	constexpr const T &z() const {
		static_assert(I <= 0, "I must satisfy I <= 0");
		return detail::iir_input_z(*this, std::integral_constant<int, I>{});
	}

private:
	template <class U, int M, int I>
	friend constexpr const U &detail::iir_input_z(const iir_input<U, M> &, std::integral_constant<int, I>);

	template <class U, int M>
	friend constexpr const U &detail::iir_input_z(const iir_input<U, M> &, std::integral_constant<int, 0>);

	T top_;
	const zsequence<T, N> &seq_;
};

template <class T>
struct iir_input<T, 0> {

	constexpr explicit iir_input(T top)
		: top_(top)
	{}

	template <int I>
	constexpr const T &z() const {
		static_assert(I == 0, "I must equal to zero");
		return top_;
	}

private:
	T top_;
};

template <class T, int N>
struct iir_output {

	constexpr explicit iir_output(const zsequence<T, N> &output)
		: seq_(output)
	{}

	template <int I>
	constexpr decltype(auto) z() const {
		static_assert(I <= 0, "I must satisfy I <= 0");
		return detail::iir_output_z(*this, std::integral_constant<int, I>{});
	}

private:
	template <class U, int M, int I>
	friend constexpr const U &detail::iir_output_z(const iir_output<U, M> &, std::integral_constant<int, I>);

	template <class U, int M>
	friend constexpr const U &detail::iir_output_z(const iir_output<U, M> &, std::integral_constant<int, 0>);

	const zsequence<T, N> &seq_;
};

template <class T>
struct iir_output<T, 0> {

	constexpr explicit iir_output() = default;

	template <int I>
	constexpr iir_result_not_set<T> z() const {
		static_assert(I == 0, "I must equal to zero");
		return iir_result_not_set<T>{};
	}
};

template <class T, int M, int N, class F>
struct iir_filter {

	constexpr explicit iir_filter(F func)
		: func_(func), input_seq_{}, output_seq_{}
	{}

	constexpr T operator()(T input) {
		iir_result_set<T> result = func_(iir_input<T, M>(input, input_seq_), iir_output<T, N>(output_seq_));

		input_seq_.push(input);
		output_seq_.push(result.value);

		return result.value;
	}

private:
	F func_;
	zsequence<T, M> input_seq_;
	zsequence<T, N> output_seq_;
};

template <class T, int N, class F>
struct iir_filter<T, 0, N, F> {

	constexpr explicit iir_filter(F func)
		: func_(func), output_seq_{}
	{}

	constexpr T operator()(T input) {
		iir_result_set<T> result = func_(iir_input<T, 0>(input), iir_output<T, N>(output_seq_));
		output_seq_.push(result.value);
		return result.value;
	}

private:
	F func_;
	zsequence<T, N> output_seq_;
};

namespace detail {
	template <class T, int M, class F>
	std::enable_if_t<
		std::is_same<
			decltype(std::declval<F>()(std::declval<iir_input<T, M>>(), std::declval<iir_output<T, 0>>())),
			iir_result_set<T>
		>::value, T
	>
	constexpr iir_filter_fir_call(F func, iir_input<T, M> input_seq) {
		return func(input_seq, iir_output<T, 0>()).value;
	}

	template <class T, int M, class F>
	std::enable_if_t<
		std::is_same<
			decltype(std::declval<F>()(std::declval<iir_input<T, M>>()))
			, T
		>::value, T
	>
	constexpr iir_filter_fir_call(F func, iir_input<T, M> input_seq)
	{
		return func(input_seq);
	}

}

template <class T, int M, class F>
struct iir_filter<T, M, 0, F> {

	constexpr explicit iir_filter(F func)
		: func_(func), input_seq_{}
	{}

	constexpr T operator()(T input) {
		T output = detail::iir_filter_fir_call<T, M>(func_, iir_input<T, M>(input, input_seq_));
		input_seq_.push(input);
		return output;
	}

private:
	F func_;
	zsequence<T, M> input_seq_;
};

template <class T, std::size_t M, std::size_t N, class F>
constexpr iir_filter<T, M, N, F> make_iir_filter(F func) {
	return iir_filter<T, M, N, F>(func);
}

template <class T, int M, int N, class U>
struct iir_filter_from_op : iir_filter<T, M, N, U> {
	template <class ...Args>
	constexpr iir_filter_from_op(Args &&...args)
		: iir_filter<T, M, N, U>(U(std::forward<Args>(args)...))
	{}
};

}

}

#endif
