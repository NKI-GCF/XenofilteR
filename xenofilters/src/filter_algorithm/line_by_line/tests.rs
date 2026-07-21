// Read-name suffix stripping (StripReadSuffix) and the LineByLine::is_new_qname
// closure it configured have been removed. Fragment-boundary detection now
// compares FragmentState::first_qname() directly; no dedicated unit tests
// remain for this behavior.
