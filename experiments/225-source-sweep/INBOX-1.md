# INBOX 1 — DISK BUDGET. Read before you fetch anything.

**The volume is at 91% — 175 GB free — and four agents are running on it.** Your mandate authorises
downloads; it does not authorise filling the disk. A disk-full event kills 224, 222 and 221 along
with you, and 222 is holding a built 12.6 GB submission package we are about to spend a slot on.

Hard rules for Phase 2:

1. **Budget: 60 GB total for everything you fetch.** Track it and report actual usage in `RESULTS.md`.
2. **Check free space before every fetch** and abort the fetch if it would leave under **100 GB**:
   `df -k . | awk 'NR==2{print $4*1024}'`. Never start a download you have not size-checked first
   (HTTP `HEAD`, or the accession's published file size).
3. **Fetch smallest-per-panel-overlap first.** A 4 GB source covering 180 panel targets outranks a
   90 GB atlas covering 200. Your Phase 1 census should already rank them this way — follow it.
4. **Delete each raw archive as soon as you have extracted the panel-restricted table you need.**
   Keep the derived table (small), not the source (huge). Record the URL and sha256 so it is
   refetchable — provenance replaces the bytes.
5. **If the census says the decisive source is too big for the budget, say so and stop.** Name it,
   size it, and put it in `RESULTS.md` as a costed recommendation. **A named, sized, unfetched
   source is a good outcome** — Ahmet has AWS credits and we can rent a box with a big volume for
   exactly that, but that is a decision to take deliberately, not by filling this laptop at 4am.

Nothing else in the mandate changes.
