/* Where the optional account lives.
 *
 * These two values are meant to be public: the key only lets a request act as an
 * anonymous visitor, and which rows anyone may read or write is decided by the
 * database's own policies, not by this file. Leaving them empty turns the account off
 * and the application keeps everything on the machine it runs on.
 */
window.GC_CLOUD = {
  url: 'https://vgayfehlczvqxqlwibde.supabase.co',
  anonKey: 'sb_publishable_e4wuMsGoJ14FojfDE40cYQ_zRPI4HZ1'
};
