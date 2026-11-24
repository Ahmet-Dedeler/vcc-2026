# Model names (`vcc submit -m`)

Decided 2026-09-26. New arms take the **next unused name from the queue below, in order**. Do not
invent names on the fly. Get it with:

```bash
uv run python scripts/next_model_name.py            # print the next free name
uv run python scripts/next_model_name.py --claim 472a   # print it and reserve it for arm 472a
```

"Used" means: in `submissions/SUBMITTED_LEDGER.json` notes, a captured board record
(`data/board-entries/*.json`), baked into any `experiments/*/chain*.sh` / `arm*.sh`, or claimed in
`docs/model-names-claimed.tsv`. Claim the name when you write it into a chain, so two chains armed
the same day never share a name.

Theme (Ahmet's pick, 2026-09-26): quotes he saved in his Notion notebook, his own and other
people's, alternating. §14 still applies: nothing that reads as a method word (union, gate, prior,
cap, whiten, donor, stack, emission, ridge, knee, sharpen, bulk, rank, ensemble, cell-line names)
and no numbers. Append new names at the end after the same screen. Never reorder. Keep the text
exactly as written here (the helper matches it verbatim, case-insensitive).

The 09-27..10-01 chains already had their names baked in (lab-object template). The queue starts
with the first arm armed after that.

**2026-09-28 15:20Z: that exception is closed.** The lab-object names that reached the board are 09-27's and 09-28's
(the 09-28 pair only because the coordinator re-armed them verbatim after a reboot, instead of switching them). Every arm
that had not been submitted yet (457dr, 446bd2, 458d, 458g, 458, 459dr, 469dr, 461dr, 460dr, 458e) now carries a queue
name claimed at 15:18Z, in chain0929 / chain0930c / chain1001b / chain1004 / chain1005 and the re-armed 09-29 waiters.
**Re-arming an inherited waiter counts as arming it: take its name from the queue.**

## Reserved

- **Final Nov 5 entry:** `Be extremely selfish. Be extremely generous.` (not in the queue; the helper skips it).

## Queue
2026-09-29 14:10Z: #57–#85 appended from the same notebook page (its middle section; same screen, also excluding
"state", "identity", "context", "random"). **Keep ≥ 20 free names for the final window (Oct 22 – Nov 4, ~16 arms):**
if the next free number before Oct 22 would leave fewer than 20, stop and extend the queue first (531 finding 19).

2026-09-29 01:30Z: #24–#56 appended from the same Notion notebook page (`Notebook`, 118d14cd…), others' and Ahmet's alternating where the page allowed, after the same screen (no digits, no method-like words such as weight / fit / norm / regress, nothing personal, religious, political or unkind). Text trimmed to the quotable sentence, apostrophes straightened.


1. Be the change you wanna see in the world.
2. You can also commit injustice by doing nothing. (Marcus Aurelius)
3. taking action is the best action—anytime. dont let it stay in your mind, go grind!!
4. If you don't believe it or don't get it, I don't have the time to try to convince you, sorry. (Satoshi)
5. You're not respawning. What are you waiting for?
6. Even if I knew the world would end tomorrow, I'd still plant an apple tree.
7. Perhaps, we don't need new ideas, but the same ideas to be known by everyone.
8. The world is bottlenecked by the ability of organizations to turn ideas into software. (John Collison)
9. Don't compare the outputs—compare the inputs.
10. The task is not so much to see what no one has yet seen, but to think what nobody has yet thought about that which everybody sees.
11. the way you do one thing is the way you do everything
12. Don't be afraid, you can do it. (Steve Jobs)
13. Ideally, your work ethic should be so far beyond the norm that it borders on the divine—understood only by you, God, and a chosen few.
14. The future is here; it is just not evenly distributed.
15. When something truly matters and you possess the ability to create change, it becomes your moral duty to do so.
16. Microsoft has had clear competitors in the past. It's a good thing we have museums to document that. (Bill Gates)
17. Hate the chains that hold you back. Hate the doubts that whisper lies. Hate the part of you that still seeks comfort in mediocrity. Use it. Let it drive you. Become the force that reshapes the world, or be forgotten by it.
18. I would like to own all the land in this valley, as far as I can see.
19. The most dangerous worldviews are the worldviews of those who have never viewed the world.
20. I wish to see my children eating off golden plates.
21. You're behind—very behind.
22. Sun Tzu, The Art of Being A Decent Human Being
23. I see problems and I solve them.
24. The reasonable man adapts himself to the world: the unreasonable one persists in trying to adapt the world to himself. Therefore all progress depends on the unreasonable man.
25. Stop calculating, start giving.
26. I needed to replace the joy of being proven right with the joy of learning what's true.
27. Just get yourself through the door, you'll be fine afterwards.
28. There's nothing sadder than living your whole life by the ocean but never daring to set sail.
29. Great minds listen to ideas. Average minds listen to same meaningless gibberish.
30. The loser has more in common with the winner than with the person sitting on the sidelines. Both had the courage to try.
31. Being a little delusional, perhaps.
32. From those to whom much is given, much is expected.
33. Have the confidence to be confident, and I am confident that you'll succeed!
34. Nothing in the world can take the place of persistence.
35. I work hard because of everything I've done. I've done everything I've done because I work hard.
36. When you see an opportunity, just shut up and take advantage of it.
37. I am persistent but open-minded.
38. In any situation, it's never as bad as it seems; it's never as good as it seems.
39. Whoever came up with April Fools must have been a great human being.
40. Ideas are easy, execution is everything.
41. Accelerate. Accelerate. Accelerate. Don't worry, you're not as far behind as you think.
42. You do not have the luxury of failure.
43. The past exists so the future can be better.
44. Doing great thing is exhausting, but it is also exhilarating.
45. social media isnt social at this point, its just media.
46. No one cares until everyone cares.
47. LACK OF COURAGE IN THE MODERN SOFT WORLD IS A VERY BIG OPPORTUNITY.
48. It is not sufficient that I succeed. All others must succeed too.
49. Do not have a life of quiet dissatisfaction.
50. If you fall in love with the process nothing else will matter.
51. Be different. It's the only hope you have to win in this infinitely competitive world.
52. There is no excuse for ignorance.
53. People who work at grocery stores are heroes.
54. everyone wants to be the main character while no one wants to be the right hand man of the main character.
55. In the future, a moat won't be quality. It will be price.
56. im doing things many dreamt of in the past, but few now notice
57. In theory, there is no difference between theory and practice; in practice, there is.
58. Speak. Speak loudly. Life is too short to be average.
59. It's what you learn after you know it all that counts.
60. There is no rule for anything in life.
61. When a man learns to feel love, he must also bear the risk of feeling hate.
62. Situations may change, and they will. But human nature doesn't change.
63. You will regret your leniency far more than your ruthlessness.
64. Fulfillment should be the thing one maximizes for. Being happy now might not be fulfilling.
65. i hope you win all the wars you tell no one about
66. I'd rather know that I have nothing, compared to thinking I have something and having nothing.
67. you call it investing but you should call it philantrophy — after all you only give money away
68. Learning from actions is probably better than just learning directly.
69. If you always want to be saying the right thing, say the honest thing.
70. What's something others aren't willing to do, that you are?
71. You've taught me what not to be, and for that, I'm grateful.
72. I am not interested in making niche products. I want to remake entire industries.
73. The chosen ones just chose themselves.. that's the secret.. It's that simple.
74. an idea that has clarity is much better than a complex idea.
75. No matter your age you'll always wish you started younger, but today is the youngest you'll ever be again.
76. HOW BADLY DO YOU WANT TO WIN???
77. One should have a deep passion and respect and be on the pursuit of truth.
78. Creativity is one of the most beautiful things in this world.
79. We're not so different.
80. I'll certainly try my best. I hope you will too.
81. if u get anything exclusive, u must be inclusive.
82. There is a gap between who you are now and who you want to be.
83. Abundance could be the key to making sure everyone is safe.
84. Beethoven walking at streets at night, singing loudly - that is how I imagine someone genuine.
85. The stars are aligning.
