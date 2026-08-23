#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cwchar>

#if !defined(_WIN32)
#define __declspec(argument)
#define _stdcall
#define FAR
#define pascal

using BSTR = wchar_t*;
using HDC = void*;
using HGDIOBJ = void*;
using COLORREF = unsigned long;

inline constexpr int DC_BRUSH = 18;
inline constexpr int DC_PEN = 19;

inline HGDIOBJ GetStockObject(int) { return nullptr; }
inline HGDIOBJ SelectObject(HDC, HGDIOBJ object) { return object; }
inline COLORREF SetDCBrushColor(HDC, COLORREF color) { return color; }
inline COLORREF SetDCPenColor(HDC, COLORREF color) { return color; }
inline int Rectangle(HDC, int, int, int, int) { return 1; }
inline int MoveToEx(HDC, int, int, void*) { return 1; }
inline int LineTo(HDC, int, int) { return 1; }
inline int TextOutA(HDC, int, int, const char*, int) { return 1; }
inline int TextOutW(HDC, int, int, const wchar_t*, int) { return 1; }
#define TextOut TextOutW

inline unsigned int SysStringLen(const BSTR value) {
    return value == nullptr ? 0U : static_cast<unsigned int>(std::wcslen(value));
}

inline int SysReAllocString(BSTR* destination, const wchar_t* source) {
    if (destination == nullptr) {
        return 0;
    }
    const std::size_t length = source == nullptr ? 0U : std::wcslen(source);
    auto* replacement = static_cast<wchar_t*>(std::calloc(length + 1U, sizeof(wchar_t)));
    if (replacement == nullptr) {
        return 0;
    }
    if (source != nullptr) {
        std::wmemcpy(replacement, source, length);
    }
    std::free(*destination);
    *destination = replacement;
    return 1;
}
#endif
