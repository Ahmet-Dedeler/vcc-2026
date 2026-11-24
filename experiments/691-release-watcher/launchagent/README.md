# 691 LaunchAgent (DRAFT, not installed)

Purpose: survive a reboot between 10-18 18:00Z and 10-27 00:00Z without any session noticing. After login, and every
30 min, `relaunch691.sh` re-launches the release watcher (691) and the R0 adapter (692) if they were killed mid-run (no copy
running, and the log's last line is not END). Outside the window it exits at once. A run that ended on purpose is never
restarted. A login is still needed after a reboot (FileVault), and launchd starts user agents only then.

**Install only after Ahmet's explicit yes**, relayed in his own words by the coordinator. It is persistent configuration on
his Mac.

Install (10-21, after the final rehearsal):
```bash
plutil -lint experiments/691-release-watcher/launchagent/com.ahmet.vcc.watch691.plist
cp experiments/691-release-watcher/launchagent/com.ahmet.vcc.watch691.plist ~/Library/LaunchAgents/
launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.ahmet.vcc.watch691.plist
launchctl print gui/$(id -u)/com.ahmet.vcc.watch691 | head -20
```

Uninstall (once R0 is submitted, or 10-27 00:00Z at the latest; also in OCT22-CHECKLIST):
```bash
launchctl bootout gui/$(id -u)/com.ahmet.vcc.watch691
rm ~/Library/LaunchAgents/com.ahmet.vcc.watch691.plist
```

Dry test (no launchd, nothing launched): `DRY691=1 NOW691=202610220300 bash experiments/691-release-watcher/launchagent/relaunch691.sh`.
