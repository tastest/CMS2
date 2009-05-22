#ifndef DataFormats_h
#define DataFormats_h

namespace edm
{
    struct BranchKey {};
    struct EventAuxiliary {};
    struct FileFormatVersion {};
    struct LuminosityBlockAuxiliary {};
    struct ParameterSetBlob {};
    struct RunAuxiliary {};

    class BranchChildren {};
    class BranchDescription {
        public:
            struct Transients {};
    };
    class BranchID {};
    class EventEntryDescription {
        public:
            struct Transients {};
    };
    class EventEntryInfo {
        public:
            struct Transients {};
    };
    class EventID {};
    class FileID {};
    class FileIndex {
        public:
            class Element {};
            struct Transients {};
    };
    template <int I> class Hash {};
    class History {};
    class LuminosityBlockID {};
    class ModuleDescription {};
    class ProcessConfiguration {};
    class ProcessHistory {
        public:
            struct Transients {};
    };
    class ProductID {};
    class ProductRegistry {
        public:
            struct Transients {};
    };
    class RunID {};
    class Timestamp {};
    template <typename T> class Transient {};
}

#endif
